"""Paired transient wave-celerity criterion."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from hydroturing.criteria.base import (
    FAIL,
    PASS,
    CriterionIncompatibleError,
    CriterionResult,
    criterion,
    make_window,
)
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec


@dataclass(frozen=True)
class _Response:
    label: str
    centroid_s: float
    baseline_discharge_m3s: float
    baseline_depth_m: float
    c_kin_m_s: float
    peak_response_m3s: float
    baseline_cv: float


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
    """Return dQ/dA from the declared rectangular Manning rating."""
    area = width * depth
    perimeter = width + 2.0 * depth
    radius = area / perimeter
    manning_q = area * radius ** (2.0 / 3.0) * np.sqrt(slope) / roughness
    return manning_q * (
        5.0 / (3.0 * area) - 4.0 / (3.0 * width * perimeter)
    )


def _measure_state(
    run: RunResult,
    probe: ProbeSpec,
    params: dict,
    variant: str,
    state: str,
) -> _Response:
    window = make_window(run, probe)
    if "dis" not in window.table.columns or "stage" not in window.table.columns:
        return _Response(f"{variant}/{state}", *(float("nan"),) * 6)

    pulse_column = str(params.get("pulse_column", "_pulse"))
    state_column = str(params.get("state_column", "_state"))
    for name in (pulse_column, state_column):
        if name not in window.forcing.columns:
            raise ValueError(f"wave-celerity forcing is missing hidden column {name!r}")

    pulse = pd.to_numeric(window.forcing[pulse_column], errors="raise").to_numpy(float)
    labels = window.forcing[state_column].astype(str).to_numpy()
    discharge = pd.to_numeric(window.table["dis"], errors="coerce").to_numpy(float)
    stage = pd.to_numeric(window.table["stage"], errors="coerce").to_numpy(float)
    if not np.isfinite(pulse).all():
        raise ValueError("wave-celerity pulse annotation must be finite")
    if not np.isfinite(discharge).all() or not np.isfinite(stage).all():
        return _Response(f"{variant}/{state}", *(float("nan"),) * 6)

    block = np.flatnonzero(labels == state)
    if len(block) == 0 or np.any(np.diff(block) != 1):
        raise ValueError(f"state {state!r} must occupy one contiguous scored block")
    marked = block[pulse[block] > 0.0]
    if len(marked) == 0 or np.any(np.diff(marked) != 1):
        raise ValueError(f"state {state!r} needs one contiguous positive perturbation")
    first = int(marked[0])

    baseline_hours = _positive(params, "baseline_hours", 6.0)
    baseline_steps = max(2, int(round(baseline_hours / (24.0 * window.dt_days))))
    if first - int(block[0]) < baseline_steps:
        raise ValueError(
            f"state {state!r} has only {first - int(block[0])} pre-pulse rows; "
            f"{baseline_steps} are required"
        )

    baseline_q = discharge[first - baseline_steps:first]
    baseline_stage = stage[first - baseline_steps:first]
    q0 = float(np.median(baseline_q))
    stage0 = float(np.median(baseline_stage))
    q_scale = max(abs(q0), np.finfo(float).eps)
    baseline_cv = float(np.std(baseline_q) / q_scale)
    max_baseline_cv = _positive(params, "max_baseline_cv", 0.005, allow_zero=True)
    if baseline_cv > max_baseline_cv:
        raise CriterionIncompatibleError(
            f"{variant}/{state} did not reach a steady hydraulic baseline: "
            f"discharge CV {baseline_cv:.3%} > {max_baseline_cv:.3%}"
        )

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
        raise ValueError("wave_celerity_bounds requires a rectangular cross-section")

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
        return _Response(
            f"{variant}/{state}",
            float("nan"),
            q0,
            depth,
            float("nan"),
            float("nan"),
            baseline_cv,
        )

    c_kin = _rectangular_celerity(depth, width, slope, roughness)
    end = int(block[-1]) + 1
    response = np.clip(discharge[first:end] - q0, 0.0, None)
    peak = float(np.max(response))
    min_peak_fraction = _positive(params, "min_peak_fraction", 0.005, allow_zero=True)
    if peak < min_peak_fraction * q0:
        return _Response(
            f"{variant}/{state}",
            float("nan"),
            q0,
            depth,
            c_kin,
            peak,
            baseline_cv,
        )

    floor_fraction = _positive(
        params, "response_floor_fraction", 1.0e-4, allow_zero=True
    )
    response[response < floor_fraction * peak] = 0.0
    total = float(response.sum())
    if total <= 0.0:
        return _Response(
            f"{variant}/{state}",
            float("nan"),
            q0,
            depth,
            c_kin,
            peak,
            baseline_cv,
        )

    centres_s = (np.arange(first, end, dtype=float) + 0.5) * window.dt_days * 86400.0
    centroid_s = float(np.dot(centres_s, response) / total)
    return _Response(
        f"{variant}/{state}",
        centroid_s,
        q0,
        depth,
        c_kin,
        peak,
        baseline_cv,
    )


def _same_reach_except_length(short: RunResult, long: RunResult) -> None:
    keys = ("width_m", "cross_section_shape", "bed_elevation_m", "slope", "manning_n")
    mismatched = [key for key in keys if short.case.static.get(key) != long.case.static.get(key)]
    if mismatched:
        raise ValueError(
            "short/long variants must share reach geometry except length; "
            f"mismatch in {mismatched}"
        )


@criterion("wave_celerity_bounds", paired=True)
def wave_celerity_bounds(
    runs: dict[str, RunResult],
    probe: ProbeSpec,
    params: dict,
) -> CriterionResult:
    """Require positive, Manning-consistent and state-responsive wave speed."""
    expected_variants = ("short", "long")
    if tuple(probe.variants) != expected_variants or set(runs) != set(expected_variants):
        raise ValueError("wave-celerity probe requires exactly short and long variants")
    short_run, long_run = runs["short"], runs["long"]
    _same_reach_except_length(short_run, long_run)

    short_length = float(short_run.case.static["reach_length_m"])
    long_length = float(long_run.case.static["reach_length_m"])
    if not (
        np.isfinite(short_length)
        and np.isfinite(long_length)
        and 0.0 < short_length < long_length
    ):
        raise ValueError("wave-celerity requires finite positive short < long lengths")
    delta_x = 0.5 * (long_length - short_length)

    tolerance = _positive(params, "tolerance", 0.05, allow_zero=True)
    pair_tolerance = _positive(params, "baseline_pair_tolerance", 0.01, allow_zero=True)
    ordering_margin = _positive(params, "ordering_margin", 0.02, allow_zero=True)
    states = [str(value) for value in params.get("states", ["low", "medium", "high"])]

    failures: list[str] = []
    diagnostics: dict[str, Any] = {"states": {}}
    observed: list[float] = []

    for state in states:
        short = _measure_state(short_run, probe, params, "short", state)
        long = _measure_state(long_run, probe, params, "long", state)
        values = [
            short.centroid_s,
            short.baseline_discharge_m3s,
            short.baseline_depth_m,
            short.c_kin_m_s,
            short.peak_response_m3s,
            long.centroid_s,
            long.baseline_discharge_m3s,
            long.baseline_depth_m,
            long.c_kin_m_s,
            long.peak_response_m3s,
        ]
        if not np.isfinite(values).all():
            failures.append(f"{state} has no finite measurable transient response")
            observed.append(float("nan"))
            continue

        q_scale = max(
            abs(short.baseline_discharge_m3s),
            abs(long.baseline_discharge_m3s),
            np.finfo(float).eps,
        )
        pair_q_mismatch = abs(
            short.baseline_discharge_m3s - long.baseline_discharge_m3s
        ) / q_scale
        if pair_q_mismatch > pair_tolerance:
            raise CriterionIncompatibleError(
                f"{state} short/long base discharge differs by "
                f"{pair_q_mismatch:.2%} > {pair_tolerance:.2%}; "
                "the common generation delay did not cancel"
            )

        delta_t = long.centroid_s - short.centroid_s
        if delta_t <= 0.0:
            c_obs = float("nan")
            residual = float("inf")
            failures.append(
                f"{state} does not propagate downstream: "
                f"long-short centroid difference {delta_t / 3600.0:.3f} h"
            )
        else:
            c_obs = delta_x / delta_t
            c_kin = 0.5 * (short.c_kin_m_s + long.c_kin_m_s)
            residual = abs(c_obs - c_kin) / c_kin
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
            "c_kin_m_s": 0.5 * (short.c_kin_m_s + long.c_kin_m_s),
            "relative_residual": residual,
            "short_baseline_discharge_m3s": short.baseline_discharge_m3s,
            "long_baseline_discharge_m3s": long.baseline_discharge_m3s,
            "short_baseline_depth_m": short.baseline_depth_m,
            "long_baseline_depth_m": long.baseline_depth_m,
            "pair_discharge_mismatch": pair_q_mismatch,
            "short_baseline_cv": short.baseline_cv,
            "long_baseline_cv": long.baseline_cv,
        }

    if len(observed) == len(states) and np.isfinite(observed).all():
        for left_state, right_state, left, right in zip(
            states[:-1], states[1:], observed[:-1], observed[1:]
        ):
            required_gap = ordering_margin * max(abs(left), abs(right))
            if right - left <= required_gap:
                failures.append(
                    f"celerity does not increase materially from {left_state} to "
                    f"{right_state}: {left:.3f} -> {right:.3f} m/s "
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
