"""Transient flood-wave celerity against a Manning/kinematic expectation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from hydroturing.criteria.base import FAIL, PASS, CriterionResult, criterion, make_window
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec


@dataclass(frozen=True)
class _Measurement:
    state: str
    c_obs: float
    c_kin: float
    residual: float
    dt_hours: float
    short_centroid_h: float
    long_centroid_h: float


def _positive(static: dict[str, Any], key: str, variant: str) -> float:
    try:
        value = float(static[key])
    except (KeyError, TypeError, ValueError):
        raise ValueError(f"variant '{variant}' needs numeric static '{key}'") from None
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"variant '{variant}' needs positive finite static '{key}'")
    return value


def _normal_depth(q: float, width: float, slope: float, manning_n: float) -> float:
    """Invert rectangular Manning discharge by monotone bisection."""
    if q <= 0:
        raise ValueError("wave celerity needs positive base discharge")

    def discharge(h: float) -> float:
        area = width * h
        perimeter = width + 2.0 * h
        radius = area / perimeter
        return (1.0 / manning_n) * area * radius ** (2.0 / 3.0) * slope ** 0.5

    lo, hi = 1.0e-6, 1.0
    while discharge(hi) < q:
        hi *= 2.0
        if hi > 1.0e4:
            raise ValueError("could not bracket Manning normal depth")
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if discharge(mid) < q:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _kinematic_celerity(q: float, static: dict[str, Any], variant: str) -> float:
    shape = str(static.get("cross_section_shape", "")).strip().lower()
    if shape != "rectangular":
        raise ValueError(
            f"variant '{variant}' needs cross_section_shape='rectangular'"
        )
    width = _positive(static, "width_m", variant)
    slope = _positive(static, "slope", variant)
    manning_n = _positive(static, "manning_n", variant)
    h = _normal_depth(q, width, slope, manning_n)
    area = width * h
    perimeter = width + 2.0 * h
    # Q(A) = n^-1 A (A/P)^(2/3) S^(1/2).
    # Hence dQ/dA = Q [5/(3A) - 4/(3 b P)] for a rectangular section.
    celerity = q * (5.0 / (3.0 * area) - 4.0 / (3.0 * width * perimeter))
    if not np.isfinite(celerity) or celerity <= 0:
        raise ValueError(f"variant '{variant}' produced invalid kinematic celerity")
    return float(celerity)


def _centroid(run: RunResult, probe: ProbeSpec, params: dict, variant: str) -> tuple[float, float]:
    window = make_window(run, probe)
    if "dis" not in window.table:
        raise ValueError(f"variant '{variant}' does not report dis")
    discharge = pd.to_numeric(window.table["dis"], errors="coerce").to_numpy(float)
    if len(discharge) != len(window.forcing) or not np.isfinite(discharge).all():
        raise ValueError(f"variant '{variant}' has invalid discharge output")

    marker = str(params.get("event_column", "_pulse"))
    if marker not in window.forcing:
        raise ValueError(f"wave celerity needs hidden forcing column '{marker}'")
    event = pd.to_numeric(window.forcing[marker], errors="coerce").to_numpy(float)
    indices = np.flatnonzero(event > 0)
    if not len(indices):
        raise ValueError("wave-celerity scored window has no pulse")
    first, last = int(indices[0]), int(indices[-1])
    if np.any(np.diff(indices) != 1):
        raise ValueError("wave-celerity pulse must be contiguous")

    baseline_hours = float(params.get("baseline_hours", 12.0))
    baseline_steps = max(2, int(round(baseline_hours / (window.dt_days * 24.0))))
    if first < baseline_steps:
        raise ValueError("wave-celerity pulse has too little pre-event baseline")
    base = float(np.mean(discharge[first - baseline_steps:first]))
    response = np.clip(discharge - base, 0.0, None)
    # Ignore numerical settling before the perturbation.
    response[:first] = 0.0
    total = float(response.sum())
    min_fraction = float(params.get("min_response_fraction", 1.0e-4))
    pulse_depth = float(np.sum(event) * window.dt_days)
    if total <= 0 or pulse_depth <= 0:
        raise ValueError(f"variant '{variant}' has no measurable positive response")
    # Use response magnitude relative to the base flow as a simple non-degeneracy guard.
    if float(response.max()) / max(abs(base), 1.0e-12) < min_fraction:
        raise ValueError(f"variant '{variant}' response is too small to time")

    centres_h = (np.arange(len(response), dtype=float) + 0.5) * window.dt_days * 24.0
    centroid_h = float(np.dot(centres_h, response) / total)
    return centroid_h, base


def _pair(runs: dict[str, RunResult], probe: ProbeSpec, params: dict, state: str) -> _Measurement:
    short_name, long_name = f"{state}_short", f"{state}_long"
    if short_name not in runs or long_name not in runs:
        raise ValueError(f"wave celerity needs variants {short_name} and {long_name}")
    short, long = runs[short_name], runs[long_name]

    # Within a hydraulic state, forcing and every static property except length
    # must be byte-identical.
    if not short.case.forcing.reset_index(drop=True).equals(
        long.case.forcing.reset_index(drop=True)
    ):
        raise ValueError(f"{state} short/long variants must have identical forcing")
    keys = set(short.case.static) | set(long.case.static)
    for key in keys - {"reach_length_m"}:
        if short.case.static.get(key) != long.case.static.get(key):
            raise ValueError(f"{state} short/long variants differ in static '{key}'")

    length_s = _positive(short.case.static, "reach_length_m", short_name)
    length_l = _positive(long.case.static, "reach_length_m", long_name)
    if length_l <= length_s:
        raise ValueError(f"{state} long reach must exceed short reach")

    t_short, q_short = _centroid(short, probe, params, short_name)
    t_long, q_long = _centroid(long, probe, params, long_name)
    if not np.isclose(q_short, q_long, rtol=0.01, atol=1.0e-9):
        raise ValueError(f"{state} short/long base discharges disagree")
    dt_hours = t_long - t_short
    timing_floor = float(params.get("timing_floor_hours", 0.05))
    if not np.isfinite(dt_hours) or dt_hours <= timing_floor:
        return _Measurement(state, -1.0, _kinematic_celerity(q_short, short.case.static, short_name),
                            float("inf"), dt_hours, t_short, t_long)

    dx = 0.5 * (length_l - length_s)
    c_obs = dx / (dt_hours * 3600.0)
    c_kin = _kinematic_celerity(q_short, short.case.static, short_name)
    residual = (c_obs - c_kin) / c_kin
    return _Measurement(state, float(c_obs), c_kin, float(residual),
                        float(dt_hours), t_short, t_long)


@criterion("wave_celerity_bounds", paired=True)
def wave_celerity_bounds(
    runs: dict[str, RunResult], probe: ProbeSpec, params: dict
) -> CriterionResult:
    """Celerity is downstream, near dQ/dA, and rises with hydraulic state."""
    states = tuple(params.get("states", ["low", "medium", "high"]))
    if len(states) < 2:
        raise ValueError("wave_celerity_bounds needs at least two hydraulic states")
    expected = {f"{s}_{side}" for s in states for side in ("short", "long")}
    if set(runs) != expected:
        raise ValueError(
            f"wave-celerity variants mismatch: expected {sorted(expected)}, got {sorted(runs)}"
        )

    tolerance = float(params.get("relative_tolerance", 0.05))
    separation = float(params.get("ordering_margin_fraction", 0.0))
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("relative_tolerance must be finite and non-negative")
    if not np.isfinite(separation) or separation < 0:
        raise ValueError("ordering_margin_fraction must be finite and non-negative")

    measured = [_pair(runs, probe, params, state) for state in states]
    failures: list[str] = []
    worst = 0.0
    diagnostics: dict[str, Any] = {}
    for item in measured:
        if item.c_obs <= 0:
            failures.append(f"{item.state}: propagation is not downstream/resolved")
            violation = 1.0
        else:
            violation = max(0.0, abs(item.residual) - tolerance)
            if abs(item.residual) > tolerance + 1.0e-12:
                failures.append(
                    f"{item.state}: c_obs={item.c_obs:.3f} m/s differs from "
                    f"c_kin={item.c_kin:.3f} m/s by {100*item.residual:+.1f}%"
                )
        worst = max(worst, violation)
        diagnostics[item.state] = {
            "c_obs_m_s": item.c_obs,
            "c_kin_m_s": item.c_kin,
            "relative_residual": item.residual,
            "short_centroid_h": item.short_centroid_h,
            "long_centroid_h": item.long_centroid_h,
            "delta_t_h": item.dt_hours,
        }

    c = [m.c_obs for m in measured]
    ordering_failures = []
    for left, right, left_name, right_name in zip(c, c[1:], states, states[1:]):
        required = max(0.0, separation * max(abs(left), abs(right)))
        if right - left <= required:
            ordering_failures.append(
                f"{left_name}->{right_name}: {left:.3f}->{right:.3f} m/s "
                f"does not exceed the {required:.3f} m/s resolution margin"
            )
    if ordering_failures:
        failures.extend(ordering_failures)
        worst = max(worst, 1.0)

    message = (
        "celerity is downstream, within the Manning/kinematic allowance, and "
        + " < ".join(f"{s}={m.c_obs:.3f} m/s" for s, m in zip(states, measured))
        if not failures else "; ".join(failures)
    )
    return CriterionResult(
        name="wave_celerity_bounds",
        status=PASS if not failures else FAIL,
        value=float(worst),
        threshold=0.0,
        message=message,
        diagnostics={
            "relative_tolerance": tolerance,
            "ordering_margin_fraction": separation,
            "states": diagnostics,
        },
    )
