"""Regression tests for momentum/wave-celerity-bounds."""

from __future__ import annotations

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
STATES = ("low", "medium", "high")
STATE_Q = (8.0, 16.0, 32.0)
EVENT_COLUMNS = {
    "low": "_pulse_low",
    "medium": "_pulse_medium",
    "high": "_pulse_high",
}


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
            "width_m", "cross_section_shape", "slope", "manning_n",
            "reach_length_m",
        ),
        generator="generate.py",
        n_seeds=1,
        timestep="PT1H",
        period_years=24.0 / 365.0,
        period_days=24.0,
        spinup_days=2,
        max_output_mb=2.0,
        max_runtime_s=300.0,
        variants=("short", "long"),
        criteria=(Criterion("wave_celerity_bounds", {}),),
        must_pass=("reference_saint_venant",),
        must_fail={"reference_fixed_celerity": "wave_celerity_bounds"},
        provenance="synthetic",
        path=ROOT / "probes" / "momentum" / "wave-celerity-bounds",
    )


def _params(**overrides) -> dict:
    params = {
        "states": list(STATES),
        "event_columns": dict(EVENT_COLUMNS),
        "relative_tolerance": 0.05,
        "baseline_hours": 12,
        "response_hours": 72,
        "timing_floor_hours": 0.05,
        "ordering_margin_fraction": 0.0,
    }
    params.update(overrides)
    return params


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


def _c_kin(q: float, width=75.0, slope=0.0012, n=0.033) -> float:
    h = _normal_depth(q, width, slope, n)
    a = width * h
    p = width + 2*h
    return q * (5/(3*a) - 4/(3*width*p))


def _synthetic_runs(celerities: tuple[float, float, float] | None = None) -> dict[str, RunResult]:
    if celerities is None:
        celerities = tuple(_c_kin(q) for q in STATE_Q)

    probe = _probe()
    spin = 48
    block = 96
    n = spin + 3 * block
    time = pd.date_range("2001-01-01", periods=n, freq="h")
    forcing = pd.DataFrame({
        "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "pr": np.zeros(n),
        "tas": np.full(n, 15.0),
        "pet": np.zeros(n),
        **{column: np.zeros(n) for column in EVENT_COLUMNS.values()},
    })

    base_by_row = np.full(n, STATE_Q[0], dtype=float)
    pulse_centres: dict[str, float] = {}
    for i, (state, q) in enumerate(zip(STATES, STATE_Q)):
        start = spin + i * block
        stop = start + block
        base_by_row[start:stop] = q
        pulse_start = start + 24
        forcing.loc[pulse_start:pulse_start + 5, EVENT_COLUMNS[state]] = 1.0
        pulse_centres[state] = pulse_start + 3.0

    # The forcing itself is irrelevant to these direct criterion fixtures
    # except for the hidden event labels and short/long equality.
    forcing["pr"] = base_by_row * 0.864

    runs = {}
    for side, length in (("short", 4000.0), ("long", 20000.0)):
        dis = base_by_row.copy()
        centres = np.arange(n, dtype=float) + 0.5
        for state, q, celerity in zip(STATES, STATE_Q, celerities):
            travel_h = length / celerity / 3600.0
            response = np.exp(
                -0.5 * ((centres - (pulse_centres[state] + travel_h)) / 1.2) ** 2
            )
            dis += 0.2 * q * response

        case = Case(
            probe_id=f"{probe.id}@{side}",
            seed=1,
            forcing=forcing.copy(),
            static={
                "area_km2": 100.0,
                "width_m": 75.0,
                "cross_section_shape": "rectangular",
                "bed_elevation_m": 50.0,
                "slope": 0.0012,
                "manning_n": 0.033,
                "reach_length_m": length,
            },
            spinup_steps=spin,
            timestep="PT1H",
        )
        table = pd.DataFrame({"time": forcing["time"], "dis": dis})
        runs[side] = RunResult(case, table, {"status": "ok"}, 0.0)
    return runs


def test_criterion_is_registered_as_paired():
    assert is_paired("wave_celerity_bounds")


def test_non_rectangular_geometry_is_rejected():
    runs = _synthetic_runs()
    for side in ("short", "long"):
        run = runs[side]
        static = dict(run.case.static)
        static["cross_section_shape"] = "trapezoidal"
        runs[side] = RunResult(
            Case(
                probe_id=run.case.probe_id,
                seed=run.case.seed,
                forcing=run.case.forcing,
                static=static,
                spinup_steps=run.case.spinup_steps,
                timestep=run.case.timestep,
            ),
            run.table,
            run.meta,
            run.wall_seconds,
        )
    with pytest.raises(ValueError, match="rectangular"):
        get("wave_celerity_bounds")(runs, _probe(), _params())


def test_registered_probe_and_references_are_compatible():
    probe = registry.find_probe("momentum/wave-celerity-bounds")
    case = build_case(probe, gate_seeds(probe.id, 1)[0], "short")
    for name in ("reference_saint_venant", "reference_fixed_celerity"):
        assert compatibility_issues(registry.find_model(name), probe, case) == []



def test_generator_pairs_change_only_reach_length():
    probe = registry.find_probe("momentum/wave-celerity-bounds")
    seed = gate_seeds(probe.id, 1)[0]
    short = build_case(probe, seed, "short")
    long = build_case(probe, seed, "long")
    pd.testing.assert_frame_equal(short.forcing, long.forcing, check_exact=True)
    s, l = dict(short.static), dict(long.static)
    assert l.pop("reach_length_m") > s.pop("reach_length_m")
    assert s == l




def test_pair_rejects_any_static_change_besides_reach_length():
    runs = _synthetic_runs()
    long = runs["long"]
    static = dict(long.case.static)
    static["width_m"] = float(static["width_m"]) * 1.01
    runs["long"] = RunResult(
        Case(
            probe_id=long.case.probe_id,
            seed=long.case.seed,
            forcing=long.case.forcing,
            static=static,
            spinup_steps=long.case.spinup_steps,
            timestep=long.case.timestep,
        ),
        long.table,
        long.meta,
        long.wall_seconds,
    )
    with pytest.raises(ValueError, match="differ in static 'width_m'"):
        get("wave_celerity_bounds")(runs, _probe(), _params())


def test_pair_rejects_a_shifted_base_state():
    runs = _synthetic_runs()
    long = runs["long"]
    table = long.table.copy()
    table["dis"] = table["dis"] * 1.02
    runs["long"] = RunResult(long.case, table, long.meta, long.wall_seconds)
    with pytest.raises(ValueError, match="base discharges disagree"):
        get("wave_celerity_bounds")(runs, _probe(), _params())


def test_generator_has_three_isolated_pulses_and_response_tail():
    probe = registry.find_probe("momentum/wave-celerity-bounds")
    case = build_case(probe, gate_seeds(probe.id, 1)[0], "short")
    scored = case.after_spinup(case.forcing)
    starts = []
    for state in STATES:
        marker = EVENT_COLUMNS[state]
        indices = np.flatnonzero(scored[marker].to_numpy(float) > 0)
        assert len(indices) == 6
        assert np.all(np.diff(indices) == 1)
        starts.append(int(indices[0]))
        assert int(indices[-1]) + 72 < len(scored)
    assert np.diff(starts).tolist() == [8 * 24, 8 * 24]


def test_generator_stays_inside_wide_channel_allowance():
    probe = registry.find_probe("momentum/wave-celerity-bounds")
    for seed in gate_seeds(probe.id, probe.n_seeds):
        case = build_case(probe, seed, "short")
        width = float(case.static["width_m"])
        slope = float(case.static["slope"])
        n = float(case.static["manning_n"])
        q = 32.0
        exact = _c_kin(q, width=width, slope=slope, n=n)
        wide_depth = (q * n / (width * slope ** 0.5)) ** (3.0 / 5.0)
        wide = (5.0 / 3.0) * q / (width * wide_depth)
        assert abs(wide - exact) / exact < 0.05


def test_analytic_measurement_passes_when_celerity_matches_manning():
    result = get("wave_celerity_bounds")(
        _synthetic_runs(),
        _probe(),
        _params(),
    )
    assert result.passed, result.message


def test_zero_response_is_a_scientific_failure_not_an_exception():
    runs = _synthetic_runs()
    run = runs["long"]
    table = run.table.copy()
    flat = np.full(len(table), STATE_Q[0], dtype=float)
    flat[48:144] = STATE_Q[0]
    flat[144:240] = STATE_Q[1]
    flat[240:336] = STATE_Q[2]
    table["dis"] = flat
    runs["long"] = RunResult(run.case, table, run.meta, run.wall_seconds)
    result = get("wave_celerity_bounds")(runs, _probe(), _params())
    assert not result.passed
    assert result.diagnostics["variant"] == "long"
    assert "no measurable positive response" in result.message


def test_nonfinite_discharge_is_a_scientific_failure_not_an_exception():
    runs = _synthetic_runs()
    run = runs["short"]
    table = run.table.copy()
    table.loc[100, "dis"] = np.nan
    runs["short"] = RunResult(run.case, table, run.meta, run.wall_seconds)
    result = get("wave_celerity_bounds")(runs, _probe(), _params())
    assert not result.passed
    assert result.diagnostics == {"variant": "short", "nonfinite_count": 1}


def test_state_response_catches_fixed_celerity_when_magnitude_check_is_neutralized():
    result = get("wave_celerity_bounds")(
        _synthetic_runs((1.0, 1.0, 1.0)),
        _probe(),
        _params(relative_tolerance=10.0, ordering_margin_fraction=0.05),
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
    assert outcome.verdict == "FAIL"
    assert "wave_celerity_bounds" in outcome.failing


def test_saint_venant_daily_uniform_flow_contract_is_preserved():
    model = registry.find_model("reference_saint_venant")
    assert model.supports_timestep("PT1D")
    assert model.supports_timestep("PT1H")
