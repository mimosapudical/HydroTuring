"""Hourly paired reaches for transient wave-celerity measurement."""

from __future__ import annotations

import numpy as np
import pandas as pd

STEP_HOURS = 1
SPINUP_DAYS = 2
PERIOD_DAYS = 5
N_STEPS = (SPINUP_DAYS + PERIOD_DAYS) * 24
AREA_KM2 = 100.0

STATE_Q = {"low": 8.0, "medium": 16.0, "high": 32.0}
LENGTHS = {"short": 4000.0, "long": 12000.0}


def _effective_mm_day(q_m3s: float) -> float:
    return q_m3s * 86400.0 / (1.0e-3 * AREA_KM2 * 1.0e6)


def generate(seed: int, variant: str = "low_short") -> tuple[pd.DataFrame, dict]:
    try:
        state, length_kind = variant.split("_", 1)
        base_q = STATE_Q[state]
        reach_length = LENGTHS[length_kind]
    except (ValueError, KeyError):
        raise ValueError(
            f"unknown variant {variant!r}; expected state_length from "
            f"{tuple(STATE_Q)} x {tuple(LENGTHS)}"
        ) from None

    rng = np.random.default_rng(seed)
    width_m = float(rng.uniform(18.0, 28.0))
    slope = float(rng.uniform(8.0e-4, 2.0e-3))
    manning_n = float(rng.uniform(0.028, 0.038))
    bed = float(rng.uniform(40.0, 120.0))

    # A small six-hour pulse around a steady positive base flow.  The same
    # forcing is used for short and long variants of one state.  The hidden
    # marker is removed before the model sees forcing.csv.
    effective = np.full(N_STEPS, _effective_mm_day(base_q))
    pulse = np.zeros(N_STEPS)
    start = SPINUP_DAYS * 24 + 24
    duration = 6
    pulse[start:start + duration] = 0.05 * effective[start:start + duration]
    pr = effective + pulse
    tas = np.full(N_STEPS, 15.0)
    pet = np.zeros(N_STEPS)

    time = pd.date_range("2001-01-01", periods=N_STEPS, freq="h")
    forcing = pd.DataFrame({
        "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "pr": pr,
        "tas": tas,
        "pet": pet,
        "_pulse": pulse,
    })
    static = {
        "area_km2": AREA_KM2,
        "width_m": width_m,
        "cross_section_shape": "rectangular",
        "bed_elevation_m": bed,
        "slope": slope,
        "manning_n": manning_n,
        "reach_length_m": reach_length,
    }
    return forcing, static
