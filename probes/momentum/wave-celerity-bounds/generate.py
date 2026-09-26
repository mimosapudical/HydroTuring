"""Hourly paired reaches for transient wave-celerity measurement."""

from __future__ import annotations

import numpy as np
import pandas as pd

STEP_HOURS = 1
SPINUP_DAYS = 20
STATE_BLOCK_DAYS = 8
PERIOD_DAYS = 3 * STATE_BLOCK_DAYS
N_STEPS = (SPINUP_DAYS + PERIOD_DAYS) * 24
AREA_KM2 = 100.0

STATE_Q = {"low": 8.0, "medium": 16.0, "high": 32.0}
LENGTHS = {"short": 4000.0, "long": 20000.0}
SETTLE_DAYS = 4
PULSE_HOURS = 6
PULSE_FRACTION = 0.05


def generate(seed: int, variant: str = "short") -> tuple[pd.DataFrame, dict]:
    if variant not in LENGTHS:
        raise ValueError(
            f"unknown variant {variant!r}; expected one of {tuple(LENGTHS)}"
        )
    reach_length = LENGTHS[variant]

    rng = np.random.default_rng(seed)
    width_m = float(rng.uniform(60.0, 90.0))
    slope = float(rng.uniform(8.0e-4, 2.0e-3))
    manning_n = float(rng.uniform(0.028, 0.038))
    bed = float(rng.uniform(40.0, 120.0))

    # q_in is an upstream river-boundary discharge in m3/s.  The experiment
    # therefore addresses the reach-routing operator directly: precipitation,
    # land runoff generation and evapotranspiration are deliberately inactive.
    # Hidden pulse columns are stripped before the model sees forcing.csv.
    q_in = np.full(N_STEPS, STATE_Q["low"], dtype=float)
    pulses = {state: np.zeros(N_STEPS, dtype=float) for state in STATE_Q}

    scored_start = SPINUP_DAYS * 24
    block_steps = STATE_BLOCK_DAYS * 24
    settle_steps = SETTLE_DAYS * 24
    for block, (state, q_m3s) in enumerate(STATE_Q.items()):
        start = scored_start + block * block_steps
        stop = start + block_steps
        q_in[start:stop] = q_m3s
        pulse_start = start + settle_steps
        pulse_stop = pulse_start + PULSE_HOURS
        pulse = PULSE_FRACTION * q_m3s
        q_in[pulse_start:pulse_stop] += pulse
        pulses[state][pulse_start:pulse_stop] = pulse

    time = pd.date_range("2001-01-01", periods=N_STEPS, freq="h")
    forcing = pd.DataFrame({
        "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
        # Keep the standard meteorological columns present so adapters that
        # require the base contract can run, but they carry no hydraulic drive.
        "pr": np.zeros(N_STEPS, dtype=float),
        "tas": np.full(N_STEPS, 15.0),
        "pet": np.zeros(N_STEPS, dtype=float),
        "q_in": q_in,
        "_pulse_low": pulses["low"],
        "_pulse_medium": pulses["medium"],
        "_pulse_high": pulses["high"],
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
