"""Hourly low/medium/high transient pairs for wave-celerity bounds."""

from __future__ import annotations

import numpy as np
import pandas as pd


SPINUP_DAYS = 1
PERIOD_DAYS = 2
STEPS_PER_DAY = 24
N_STEPS = (SPINUP_DAYS + PERIOD_DAYS) * STEPS_PER_DAY
SHORT_REACH_M = 10_000.0
LONG_REACH_M = 50_000.0
PULSE_START = SPINUP_DAYS * STEPS_PER_DAY + 8
PULSE_STEPS = 2
STATE_MULTIPLIERS = {"low": 1.0, "medium": 2.0, "high": 4.0}
VARIANTS = tuple(
    f"{state}_{length}"
    for state in STATE_MULTIPLIERS
    for length in ("short", "long")
)


def generate(seed: int, variant: str = "low_short") -> tuple[pd.DataFrame, dict]:
    """Generate one state/length member of the six-case transient experiment."""
    if variant not in VARIANTS:
        raise ValueError(f"unknown variant {variant!r}; expected one of {VARIANTS}")

    state, length_name = variant.rsplit("_", 1)
    rng = np.random.default_rng(seed)

    # Keep the experiment inside the mild, subcritical, friction-dominated
    # regime in which a local kinematic-wave comparison is the intended limit.
    area_km2 = float(rng.uniform(900.0, 1100.0))
    width_m = float(rng.uniform(95.0, 125.0))
    slope = float(rng.uniform(1.0e-3, 2.0e-3))
    manning_n = float(rng.uniform(0.030, 0.040))
    bed_elevation_m = float(rng.uniform(40.0, 120.0))

    pet = float(rng.uniform(1.5, 2.0))
    tas = float(rng.uniform(14.0, 18.0))
    low_effective_mm_day = float(rng.uniform(0.75, 0.95))
    effective = low_effective_mm_day * STATE_MULTIPLIERS[state]

    pr = np.full(N_STEPS, pet + effective, dtype=float)
    pulse = np.zeros(N_STEPS, dtype=float)
    # Five percent is small enough to stay in the local-wave regime but large
    # enough to dominate floating-point and finite-volume noise.
    pulse[PULSE_START:PULSE_START + PULSE_STEPS] = 0.05 * effective
    pr += pulse

    time = pd.date_range("2001-01-01", periods=N_STEPS, freq="h")
    forcing = pd.DataFrame(
        {
            "time": time.strftime("%Y-%m-%d %H:%M:%S"),
            "pr": np.round(pr, 9),
            "tas": np.full(N_STEPS, tas),
            "pet": np.full(N_STEPS, pet),
            # Hidden annotation: stripped before a model sees the forcing.
            "_pulse": np.round(pulse, 9),
        }
    )

    static = {
        "area_km2": area_km2,
        "soil_capacity_mm": 320.0,
        "canopy_capacity_mm": 2.0,
        "degree_day_factor_mm_per_C_day": 3.2,
        "baseflow_coefficient": 0.006,
        "snow_threshold_degC": 0.0,
        "latitude_deg": 38.0,
        "width_m": width_m,
        "cross_section_shape": "rectangular",
        "bed_elevation_m": bed_elevation_m,
        "slope": slope,
        "manning_n": manning_n,
        "reach_length_m": (
            SHORT_REACH_M if length_name == "short" else LONG_REACH_M
        ),
        # Repository-owned Saint-Venant reference uses this opt-in to advance
        # the transient PDE instead of its existing pseudo-steady mode.
        "transient_wave": True,
    }
    return forcing, static
