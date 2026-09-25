"""Transient wave-celerity consistency for a rectangular open channel."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion, make_window
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec


@dataclass(frozen=True)
class _Response:
    variant: str
    centroid_s: float
    baseline_discharge_m3s: float
    baseline_depth_m: float
    c_kin_m_s: float
    peak_response_m3s: float


def _positive(params: dict, key: str, default: float, *, allow_zero: bool = False) -> float:
    value = float(params.get(key, default))
    ok = np.isfinite(value) and (value >= 0.0 if allow_zero else value > 0.0)
    if not ok:
        relation = "non-negative" if allow_zero else "positive"
        raise ValueError(f"wave-celerity parameter {key!r} must be finite and {relation}")
    return value


def _rectangular_celerity(
    depth: float,
    width: float,
    slope: float,
    roughness: float,
) -> float:
    """Return dQ/dA from the declared rectangular Manning rating.

    The expected discharge is reconstructed from geometry rather than borrowed
    from the model output:

        Q_M = A R^(2/3) S^(1/2) / n
        dQ_M/dA = Q_M [5/(3A) - 4/(3 w P)].

    That independence matters: a model cannot make an inconsistent Q/stage pair
    define its own expected celerity.
    """
    area = width * depth
    perimeter = width + 2.0 * depth
    radius = area / perimeter
    manning_q = area * radius ** (2.0 / 3.0) * np.sqrt(slope) / roughness
    return manning_q * (
        5.0 / (3.0 * area) - 4.0 / (3.0 * width * perimeter)
    )


def _variant_order(runs: dict[str, RunResult], probe: ProbeSpec) -> list[str]:
    ordered = list(probe.variants) if probe.variants else list(runs)
    missing = [name for name in ordered if name not in runs]
    extra = [name for name in runs if name not in ordered]
    if missing or extra:
        raise ValueError(
            "wave-celerity variants do not match the probe declaration "
            f"(missing={missing}, extra={extra})"
        )
    return ordered


def _measure(
    run: RunResult,
    probe: ProbeSpec,
    params: dict,
    variant: str,
) -> _Response:
    window = make_window(run, probe)
    for name in ("dis", "stage"):
        if name not in window.table.columns:
            return _Response(variant, np.nan, np.nan, np.nan, np.nan, np.nan)

    pulse_column = str(params.get("pulse_column", "_pulse"))
    if pulse_column not in window.forcing.columns:
        raise ValueError(f"wave-celerity forcing is missing hidden column {pulse_column!r}")

    pulse = pd.to_numeric(window.forcing[pulse_column], errors="raise").to_numpy(float)
    discharge = pd.to_numeric(window.table["dis"], errors="coerce").to_numpy(float)
    stage = pd.to_numeric(window.table["stage"], errors="coerce").to_numpy(float)
    if not np.isfinite(pulse).all():
        raise ValueError(f"variant {variant!r} has non-finite pulse forcing")
    if not np.isfinite(discharge).all() or not np.isfinite(stage).all():
        return _Response(variant, np.nan, np.nan, np.nan, np.nan, np.nan)

    marked = np.flatnonzero(pulse > 0.0)
    if len(marked) == 0 or np.any(np.diff(marked) != 1):
        raise ValueError("wave-celerity needs one contiguous positive perturbation")
    first = int(marked[0])

    baseline_hours = _positive(params, "baseline_hours", 6.0)
    baseline_steps = max(2, int(round(baseline_hours / (24.0 * window.dt_days))))
    if first < baseline_steps:
        raise ValueError(
            f"pulse starts after {first} rows, but {baseline_steps} baseline rows are required"
        )

    q0 = float(np.median(discharge[first - baseline_steps:first]))
    stage0 = float(np.median(stage[first - baseline_steps:first]))
    static = run.case.static
    required = (
        "width_m",
        "bed_elevation_m",
        "slope",
        "manning_n",
        "cross_section_shape",
        "reach_length_m",
    )
    absent = [name for name in required if name not in static]
    if absent:
        raise ValueError(f"case is missing static value(s): {', '.join(absent)}")
    width = float(static["width_m"])
    bed = float(static["bed_elevation_m"])
    slope = float(static["slope"])
    roughness = float(static["manning_n"])
    shape = str(static["cross_section_shape"]).strip().lower()
    if shape != "rectangular":
        raise ValueError(
            "wave_celerity_bounds requires cross_section_shape='rectangular'"
        )
    depth = stage0 - bed
    numbers = np.asarray([q0, width, depth, slope, roughness], dtype=float)
    if (
        not np.isfinite(numbers).all()
        or q0 <= 0.0
        or width <= 0.0
        or depth <= 0.0
        or slope <= 0.0
        or roughness <= 0.0
    ):
        return _Response(variant, np.nan, q0, depth, np.nan, np.nan)

    c_kin = _rectangular_celerity(depth, width, slope, roughness)
    if not np.isfinite(c_kin) or c_kin <= 0.0:
        return _Response(variant, np.nan, q0, depth, c_kin, np.nan)

    response = np.clip(discharge - q0, 0.0, None)
    response[:first] = 0.0
    peak = float(np.max(response))
    min_peak_fraction = _positive(params, "min_peak_fraction", 0.005, allow_zero=True)
    if peak < min_peak_fraction * q0:
        return _Response(variant, np.nan, q0, depth, c_kin, peak)

    # A tiny numerical tail can otherwise move a centroid arbitrarily far.
    floor_fraction = _positive(params, "response_floor_fraction", 1.0e-4, allow_zero=True)
    response[response < floor_fraction * peak] = 0.0
    if float(response.sum()) <= 0.0:
        return _Response(variant, np.nan, q0, depth, c_kin, peak)

    centres_s = (np.arange(len(response), dtype=float) + 0.5) * window.dt_days * 86400.0
    centroid_s = float(np.dot(centres_s, response) / response.sum())
    return _Response(variant, centroid_s, q0, depth, c_kin, peak)


@criterion("wave_celerity_bounds", paired=True)
def wave_celerity_bounds(
    runs: dict[str, RunResult],
    probe: ProbeSpec,
    params: dict,
) -> CriterionResult:
    """Require positive, Manning-consistent and state-responsive wave speed."""
    ordered = _variant_order(runs, probe)
    states = [str(value) for value in params.get("states", ["low", "medium", "high"])]
    expected = [f"{state}_{length}" for state in states for length in ("short", "long")]
    if ordered != expected:
        raise ValueError(
            "wave-celerity variants must be ordered "
            + ", ".join(expected)
        )

    tolerance = _positive(params, "tolerance", 0.05, allow_zero=True)
    baseline_pair_tolerance = _positive(
        params, "baseline_pair_tolerance", 0.01, allow_zero=True
    )
    ordering_margin = _positive(params, "ordering_margin", 0.02, allow_zero=True)

    measurements = {
        name: _measure(runs[name], probe, params, name)
        for name in ordered
    }
    failures: list[str] = []
    diagnostics: dict[str, Any] = {"states": {}}
    observed: list[float] = []

    for state in states:
        short = measurements[f"{state}_short"]
        long = measurements[f"{state}_long"]
        invalid = []
        for item in (short, long):
            if not np.isfinite(
                [
                    item.centroid_s,
                    item.baseline_discharge_m3s,
                    item.baseline_depth_m,
                    item.c_kin_m_s,
                    item.peak_response_m3s,
                ]
            ).all():
                failures.append(
                    f"{item.variant} has no finite measurable transient response"
                )
                invalid.append(item.variant)

        if invalid:
            continue

        q_scale = max(
            abs(short.baseline_discharge_m3s),
            abs(long.baseline_discharge_m3s),
            np.finfo(float).eps,
        )
        pair_q_mismatch = abs(
            short.baseline_discharge_m3s - long.baseline_discharge_m3s
        ) / q_scale
        if pair_q_mismatch > baseline_pair_tolerance:
            raise ValueError(
                f"{state} short/long baseline discharge differs by {pair_q_mismatch:.2%}; "
                "only reach length may change inside a state pair"
            )

        short_length = float(runs[f"{state}_short"].case.static["reach_length_m"])
        long_length = float(runs[f"{state}_long"].case.static["reach_length_m"])
        if not (np.isfinite(short_length) and np.isfinite(long_length) and 0 < short_length < long_length):
            raise ValueError(f"{state} needs finite positive short < long reach lengths")
        delta_x = 0.5 * (long_length - short_length)
        delta_t = long.centroid_s - short.centroid_s
        if delta_t <= 0.0:
            failures.append(
                f"{state} response does not propagate downstream: "
                f"long-short centroid difference is {delta_t / 3600.0:.3f} h"
            )
            c_obs = float("nan")
        else:
            c_obs = delta_x / delta_t

        c_kin = 0.5 * (short.c_kin_m_s + long.c_kin_m_s)
        residual = abs(c_obs - c_kin) / c_kin if np.isfinite(c_obs) else float("inf")
        if residual > tolerance:
            failures.append(
                f"{state} celerity residual {residual:.2%} exceeds {tolerance:.2%}"
            )
        observed.append(c_obs)
        diagnostics["states"][state] = {
            "short_centroid_h": short.centroid_s / 3600.0,
            "long_centroid_h": long.centroid_s / 3600.0,
            "delta_t_h": delta_t / 3600.0,
            "delta_x_m": delta_x,
            "c_obs_m_s": c_obs,
            "c_kin_m_s": c_kin,
            "relative_residual": residual,
            "short_baseline_discharge_m3s": short.baseline_discharge_m3s,
            "long_baseline_discharge_m3s": long.baseline_discharge_m3s,
            "short_baseline_depth_m": short.baseline_depth_m,
            "long_baseline_depth_m": long.baseline_depth_m,
            "pair_discharge_mismatch": pair_q_mismatch,
        }

    if len(observed) == len(states) and np.isfinite(observed).all():
        for left_state, right_state, left, right in zip(
            states[:-1], states[1:], observed[:-1], observed[1:]
        ):
            required_gap = ordering_margin * max(abs(left), abs(right))
            if right - left <= required_gap:
                failures.append(
                    f"celerity does not increase materially from {left_state} to {right_state}: "
                    f"{left:.3f} -> {right:.3f} m/s "
                    f"(required gap > {required_gap:.3f} m/s)"
                )

    finite_residuals = [
        float(item["relative_residual"])
        for item in diagnostics["states"].values()
        if np.isfinite(item["relative_residual"])
    ]
    value = max(finite_residuals) if finite_residuals else float("inf")
    diagnostics["ordering_margin"] = ordering_margin
    diagnostics["observed_celerities_m_s"] = observed

    if failures:
        return CriterionResult(
            "wave_celerity_bounds",
            FAIL,
            "; ".join(failures),
            value=value,
            threshold=tolerance,
            diagnostics=diagnostics,
        )

    summary = ", ".join(
        f"{state} {diagnostics['states'][state]['c_obs_m_s']:.3f}/"
        f"{diagnostics['states'][state]['c_kin_m_s']:.3f} m/s"
        for state in states
    )
    return CriterionResult(
        "wave_celerity_bounds",
        PASS,
        f"observed/kinematic celerities agree and increase with state ({summary})",
        value=value,
        threshold=tolerance,
        diagnostics=diagnostics,
    )
