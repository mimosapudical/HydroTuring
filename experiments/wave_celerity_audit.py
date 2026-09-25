#!/usr/bin/env python3
"""Feasibility audit for a state-dependent wave-celerity probe.

This is deliberately *not* a registered HydroTuring probe.  It reuses the
merged `reference_saint_venant` finite-volume solver and asks one narrow
question before a proposal is opened:

Can a transient Saint-Venant reach distinguish state-dependent flood-wave
celerity from a fixed-celerity router that already passes the merged
`momentum/routing-lag-consistency` probe?

Method
------
For each seed we draw one mild, wide rectangular reach and two reach lengths.
At three steady base depths (0.5, 1.0, 2.0 m), we:
  1. solve the existing Saint-Venant model to steady state;
  2. add a small Gaussian upstream-discharge pulse;
  3. advance the existing solver in real time via its private `_advance`;
  4. measure the positive-response centroid at the model's center gauge;
  5. difference the centroid times between short and long reaches.

Because both variants receive the same input pulse, its generation time cancels.
The observed wave speed is

    c_obs = ((L_long - L_short) / 2) / (t_long - t_short)

The comparison value is the local kinematic-wave speed dQ/dA computed from the
exact rectangular Manning rating.  The Saint-Venant solver never calls that
rating to advance the transient, so this is an independent comparison path.

This audit is exploratory evidence only.  It does not freeze a production
criterion or tolerance.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOLVER_PATH = ROOT / "models" / "reference_saint_venant" / "ht_adapter.py"
BASE_DEPTHS_M = (0.5, 1.0, 2.0)
DEFAULT_SEEDS = 50
DEFAULT_N_CELLS = 64
PULSE_FRACTION = 0.05
PULSE_CENTER_HOURS = 4.0
PULSE_SIGMA_HOURS = 0.5


def load_solver():
    spec = importlib.util.spec_from_file_location("saint_venant_reference", SOLVER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {SOLVER_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def manning_discharge(depth_m: float, width_m: float, slope: float, manning_n: float) -> float:
    area = width_m * depth_m
    radius = area / (width_m + 2.0 * depth_m)
    return (1.0 / manning_n) * area * radius ** (2.0 / 3.0) * math.sqrt(slope)


def kinematic_celerity(depth_m: float, width_m: float, slope: float, manning_n: float) -> float:
    """Return dQ/dA for the exact rectangular Manning rating."""
    eps = max(1.0e-7, depth_m * 1.0e-5)
    dq_dh = (
        manning_discharge(depth_m + eps, width_m, slope, manning_n)
        - manning_discharge(depth_m - eps, width_m, slope, manning_n)
    ) / (2.0 * eps)
    return dq_dh / width_m


def transient_centroid(
    solver: Any,
    discharge_m3s: float,
    width_m: float,
    slope: float,
    manning_n: float,
    reach_length_m: float,
    *,
    n_cells: int,
    duration_h: float,
    pulse_fraction: float = PULSE_FRACTION,
) -> dict[str, float]:
    depth, unit_q, steady = solver.solve_steady_reach(
        discharge_m3s,
        width_m,
        slope,
        manning_n,
        reach_length_m,
        n_cells=n_cells,
    )
    center = n_cells // 2
    base_q = float(width_m * unit_q[center])
    dx_m = reach_length_m / n_cells

    time_s = 0.0
    times: list[float] = []
    gauge_q: list[float] = []
    while time_s < duration_h * 3600.0:
        pulse = pulse_fraction * math.exp(
            -0.5
            * (
                (time_s / 3600.0 - PULSE_CENTER_HOURS)
                / PULSE_SIGMA_HOURS
            )
            ** 2
        )
        upstream_q = discharge_m3s * (1.0 + pulse)
        depth, unit_q, dt_s = solver._advance(
            depth,
            unit_q,
            upstream_q / width_m,
            width_m,
            slope,
            manning_n,
            dx_m,
        )
        time_s += dt_s
        times.append(time_s)
        gauge_q.append(float(width_m * unit_q[center]))

    t = np.asarray(times, dtype=float)
    response = np.maximum(np.asarray(gauge_q, dtype=float) - base_q, 0.0)
    area = float(np.trapezoid(response, t))
    if not math.isfinite(area) or area <= 0.0:
        raise RuntimeError("transient produced no positive gauge response")
    centroid_s = float(np.trapezoid(response * t, t) / area)
    return {
        "centroid_s": centroid_s,
        "peak_s": float(t[int(np.argmax(response))]),
        "base_discharge_m3s": base_q,
        "max_response_m3s": float(response.max()),
        "steady_steps": int(steady["steps"]),
        "transient_steps": int(len(t)),
    }


def run_seed(solver: Any, seed: int, *, n_cells: int = DEFAULT_N_CELLS) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    width_m = float(rng.uniform(60.0, 120.0))
    slope = float(rng.uniform(3.0e-4, 1.2e-3))
    manning_n = float(rng.uniform(0.030, 0.045))
    short_length_m = float(rng.uniform(20_000.0, 30_000.0))
    long_length_m = float(rng.uniform(65_000.0, 85_000.0))

    levels = []
    for depth_m in BASE_DEPTHS_M:
        discharge = manning_discharge(depth_m, width_m, slope, manning_n)
        theory_c = kinematic_celerity(depth_m, width_m, slope, manning_n)
        velocity = discharge / (width_m * depth_m)
        froude = velocity / math.sqrt(solver.GRAVITY * depth_m)

        duration_h = max(
            36.0,
            PULSE_CENTER_HOURS
            + 2.0 * (0.5 * long_length_m / theory_c) / 3600.0
            + 8.0,
        )
        short = transient_centroid(
            solver,
            discharge,
            width_m,
            slope,
            manning_n,
            short_length_m,
            n_cells=n_cells,
            duration_h=duration_h,
        )
        long = transient_centroid(
            solver,
            discharge,
            width_m,
            slope,
            manning_n,
            long_length_m,
            n_cells=n_cells,
            duration_h=duration_h,
        )

        travel_distance_difference_m = 0.5 * (long_length_m - short_length_m)
        lag_s = long["centroid_s"] - short["centroid_s"]
        if lag_s <= 0.0:
            raise RuntimeError(f"seed {seed}: long reach did not lag short reach")
        observed_c = travel_distance_difference_m / lag_s

        levels.append(
            {
                "base_depth_m": depth_m,
                "base_discharge_m3s": discharge,
                "froude": froude,
                "theory_celerity_m_s": theory_c,
                "observed_celerity_m_s": observed_c,
                "observed_to_theory": observed_c / theory_c,
                "short_centroid_h": short["centroid_s"] / 3600.0,
                "long_centroid_h": long["centroid_s"] / 3600.0,
            }
        )

    return {
        "seed": seed,
        "width_m": width_m,
        "slope": slope,
        "manning_n": manning_n,
        "short_length_m": short_length_m,
        "long_length_m": long_length_m,
        "levels": levels,
    }


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    observed = np.asarray(
        [[level["observed_celerity_m_s"] for level in row["levels"]] for row in rows],
        dtype=float,
    )
    theory = np.asarray(
        [[level["theory_celerity_m_s"] for level in row["levels"]] for row in rows],
        dtype=float,
    )
    ratios = observed / theory
    monotonic = (observed[:, 0] < observed[:, 1]) & (observed[:, 1] < observed[:, 2])
    high_low = observed[:, 2] / observed[:, 0]

    return {
        "monotonic_pass_seeds": int(monotonic.sum()),
        "total_seeds": int(len(rows)),
        "observed_to_theory_ratio_min": float(ratios.min()),
        "observed_to_theory_ratio_median": float(np.median(ratios)),
        "observed_to_theory_ratio_max": float(ratios.max()),
        "high_to_low_celerity_ratio_min": float(high_low.min()),
        "high_to_low_celerity_ratio_median": float(np.median(high_low)),
        "high_to_low_celerity_ratio_max": float(high_low.max()),
        "max_froude": float(
            max(level["froude"] for row in rows for level in row["levels"])
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", type=int, default=DEFAULT_SEEDS)
    parser.add_argument("--cells", type=int, default=DEFAULT_N_CELLS)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    solver = load_solver()
    rows = [run_seed(solver, seed, n_cells=args.cells) for seed in range(args.seeds)]
    payload = {
        "design": {
            "seeds": args.seeds,
            "n_cells": args.cells,
            "pulse_fraction": PULSE_FRACTION,
            "pulse_center_hours": PULSE_CENTER_HOURS,
            "pulse_sigma_hours": PULSE_SIGMA_HOURS,
            "base_depths_m": list(BASE_DEPTHS_M),
            "fixed_celerity_control_m_s": [1.0, 1.0, 1.0],
        },
        "summary": summarize(rows),
        "per_seed": rows,
    }
    rendered = json.dumps(payload, indent=2)
    if args.output:
        args.output.write_text(rendered + "\n")
    print(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
