"""Regression tests for momentum/wave-celerity-bounds."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria import get, is_paired
from hydroturing.harness import build_case, compatibility_issues, run_probe
from hydroturing.protocol import Case, RunResult
from hydroturing.seeds import gate_seeds
from hydroturing.spec import Criterion, ProbeSpec

ROOT = Path(__file__).resolve().parents[1]


def _probe() -> ProbeSpec:
    return ProbeSpec(
        id="momentum/wave-celerity-bounds",
        title="test",
        law="momentum",
        track="synthetic",
        version=1,
        authors=({"name": "test", "affiliation": "test"},),
        citation="",
        requires_fluxes=("dis",),
        requires_states=(),
        requires_forcing=("pr",),
        requires_static=(
            "area_km2", "width_m", "cross_section_shape", "bed_elevation_m",
            "slope", "manning_n", "reach_length_m",
        ),
        generator="generate.py",
        n_seeds=1,
        timestep="PT1H",
        period_years=5.0 / 365.0,
        period_days=5.0,
        spinup_days=2,
        max_output_mb=2.0,
        max_runtime_s=300.0,
        variants=(
            "low_short", "low_long",
            "medium_short", "medium_long",
            "high_short", "high_long",
        ),
        criteria=(Criterion("wave_celerity_bounds", {}),),
        must_pass=("reference_saint_venant",),
        must_fail={"reference_fixed_celerity": "wave_celerity_bounds"},
        provenance="synthetic",
        path=ROOT / "probes" / "momentum" / "wave-celerity-bounds",
    )


def _normal_depth(q: float, width: float, slope: float, n: float) -> float:
    def f(h):
        a = width * h
        r = a / (width + 2*h)
        return a * r ** (2/3) * slope ** 0.5 / n
    lo, hi = 1e-6, 1.0
    while f(hi) < q:
        hi *= 2
    for _ in range(80):
        mid = (lo + hi) / 2
        if f(mid) < q:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def _c_kin(q: float, width=22.0, slope=0.0012, n=0.033) -> float:
    h = _normal_depth(q, width, slope, n)
    a = width * h
    p = width + 2*h
    return q * (5/(3*a) - 4/(3*width*p))


def _synthetic_runs(celerities=(1.1, 1.5, 2.0)) -> dict[str, RunResult]:
    probe = _probe()
    n = 168
    spin = 48
    start = 72
    pulse = np.zeros(n)
    pulse[start:start+6] = 1.0
    time = pd.date_range("2001-01-01", periods=n, freq="h")
    runs = {}
    for state, q, c in zip(("low", "medium", "high"), (8.0, 16.0, 32.0), celerities):
        for side, length in (("short", 4000.0), ("long", 12000.0)):
            forcing = pd.DataFrame({
                "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
                "pr": np.full(n, q * 0.864),
                "tas": np.full(n, 15.0),
                "pet": np.zeros(n),
                "_pulse": pulse,
            })
            # A narrow Gaussian response whose centroid is travel time L/(2c).
            travel_h = 0.5 * length / c / 3600.0
            centres = np.arange(n, dtype=float) + 0.5
            event_centroid = start + 3.0
            response = np.exp(-0.5 * ((centres - (event_centroid + travel_h)) / 1.2) ** 2)
            table = pd.DataFrame({"time": forcing["time"], "dis": q + 0.2 * q * response})
            case = Case(
                probe_id=f"{probe.id}@{state}_{side}",
                seed=1,
                forcing=forcing,
                static={
                    "area_km2": 100.0,
                    "width_m": 22.0,
                    "cross_section_shape": "rectangular",
                    "bed_elevation_m": 50.0,
                    "slope": 0.0012,
                    "manning_n": 0.033,
                    "reach_length_m": length,
                },
                spinup_steps=spin,
                timestep="PT1H",
            )
            runs[f"{state}_{side}"] = RunResult(case, table, {"status": "ok"}, 0.0)
    return runs


def test_criterion_is_registered_as_paired():
    assert is_paired("wave_celerity_bounds")


def test_registered_probe_and_references_are_compatible():
    probe = registry.find_probe("momentum/wave-celerity-bounds")
    case = build_case(probe, gate_seeds(probe.id, 1)[0], "low_short")
    for name in ("reference_saint_venant", "reference_fixed_celerity"):
        assert compatibility_issues(registry.find_model(name), probe, case) == []


def test_generator_pairs_change_only_reach_length():
    probe = registry.find_probe("momentum/wave-celerity-bounds")
    seed = gate_seeds(probe.id, 1)[0]
    for state in ("low", "medium", "high"):
        short = build_case(probe, seed, f"{state}_short")
        long = build_case(probe, seed, f"{state}_long")
        pd.testing.assert_frame_equal(short.forcing, long.forcing, check_exact=True)
        s, l = dict(short.static), dict(long.static)
        assert l.pop("reach_length_m") > s.pop("reach_length_m")
        assert s == l


def test_analytic_measurement_passes_when_celerity_matches_manning():
    expected = tuple(_c_kin(q) for q in (8.0, 16.0, 32.0))
    result = get("wave_celerity_bounds")(
        _synthetic_runs(expected),
        _probe(),
        {
            "event_column": "_pulse",
            "states": ["low", "medium", "high"],
            "relative_tolerance": 0.05,
            "baseline_hours": 12,
            "timing_floor_hours": 0.05,
        },
    )
    assert result.passed, result.message


def test_fixed_celerity_fails_state_response():
    result = get("wave_celerity_bounds")(
        _synthetic_runs((1.0, 1.0, 1.0)),
        _probe(),
        {
            "event_column": "_pulse",
            "states": ["low", "medium", "high"],
            # Isolate the ordering assertion in this unit test.
            "relative_tolerance": 10.0,
            "baseline_hours": 12,
            "timing_floor_hours": 0.05,
        },
    )
    assert not result.passed
    assert "low->medium" in result.message or "medium->high" in result.message


def test_registered_fixed_celerity_control_trips_gate():
    probe = registry.find_probe("momentum/wave-celerity-bounds")
    outcome = run_probe(
        registry.find_model("reference_fixed_celerity"),
        probe,
        gate_seeds(probe.id, 1),
    )
    assert outcome.verdict == "fail"
    assert "wave_celerity_bounds" in outcome.failing


def test_saint_venant_daily_uniform_flow_contract_is_preserved():
    model = registry.find_model("reference_saint_venant")
    assert model.supports_timestep("PT1D")
    assert model.supports_timestep("PT1H")
