"""Focused tests for transient wave-celerity bounds."""

from __future__ import annotations

import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria import is_paired
from hydroturing.harness import build_case, compatibility_issues, run_probe
from hydroturing.scoring import FAIL, PASS
from hydroturing.seeds import gate_seeds


PROBE_ID = "momentum/wave-celerity-bounds"
STATES = ("low", "medium", "high")
STANDARD_BASELINES = (
    "reference_bucket",
    "flex_lumped",
    "flex_topo",
    "sacsma_snow17",
)


@pytest.fixture(scope="module")
def probe():
    return registry.find_probe(PROBE_ID)


def test_criterion_is_paired_and_is_the_probe_headline(probe):
    assert is_paired("wave_celerity_bounds")
    assert probe.variants == ("short", "long")
    assert probe.headline == ("wave_celerity_bounds",)


def test_probe_declares_every_geometry_input_used_by_the_verdict(probe):
    assert probe.requires_fluxes == ("dis",)
    assert probe.requires_diagnostics == ("stage",)
    assert probe.requires_static == (
        "width_m",
        "cross_section_shape",
        "bed_elevation_m",
        "slope",
        "manning_n",
        "reach_length_m",
    )


@pytest.mark.parametrize("seed", gate_seeds(PROBE_ID, 3))
def test_short_long_runs_are_identical_except_reach_length(probe, seed):
    short = build_case(probe, seed, "short")
    long = build_case(probe, seed, "long")
    pd.testing.assert_frame_equal(short.forcing, long.forcing, check_exact=True)
    assert set(short.static) == set(long.static)
    for key in short.static:
        if key != "reach_length_m":
            assert short.static[key] == long.static[key]
    assert short.static["reach_length_m"] == pytest.approx(10_000.0)
    assert long.static["reach_length_m"] == pytest.approx(50_000.0)


@pytest.mark.parametrize("seed", gate_seeds(PROBE_ID, 3))
def test_each_state_has_settling_time_and_one_small_pulse(probe, seed):
    case = build_case(probe, seed, "short")
    scored = case.after_spinup(case.forcing)
    bases = []
    for state in STATES:
        block = scored.loc[scored["_state"] == state]
        assert len(block) == 72
        pulse = block["_pulse"] > 0.0
        assert int(pulse.sum()) == 2
        first_pulse = int(pulse.to_numpy().nonzero()[0][0])
        assert first_pulse == 36
        base = float((block["pr"] - block["pet"] - block["_pulse"]).iloc[0])
        bases.append(base)
        assert float(block["_pulse"].max()) == pytest.approx(0.05 * base)
    assert bases[1] == pytest.approx(2.0 * bases[0])
    assert bases[2] == pytest.approx(4.0 * bases[0])


def test_only_process_aware_references_are_compatible(probe):
    case = build_case(probe, gate_seeds(PROBE_ID, 1)[0], "short")
    assert compatibility_issues(
        registry.find_model("reference_saint_venant"), probe, case
    ) == []
    assert compatibility_issues(
        registry.find_model("reference_fixed_celerity"), probe, case
    ) == []
    for model_name in STANDARD_BASELINES:
        issues = compatibility_issues(registry.find_model(model_name), probe, case)
        assert any("reach_length_m" in issue for issue in issues), (model_name, issues)


def test_rectangular_celerity_is_the_analytic_manning_derivative():
    module = __import__(
        "hydroturing.criteria.wave_celerity",
        fromlist=["_rectangular_celerity"],
    )
    width = 100.0
    depth = 0.8
    slope = 0.0015
    roughness = 0.035
    exact = module._rectangular_celerity(depth, width, slope, roughness)
    area = width * depth

    def rating(a):
        d = a / width
        r = a / (width + 2.0 * d)
        return a * r ** (2.0 / 3.0) * slope ** 0.5 / roughness

    eps = 1.0e-4 * area
    finite_difference = (rating(area + eps) - rating(area - eps)) / (2.0 * eps)
    assert exact == pytest.approx(finite_difference, rel=1.0e-8)


def test_fixed_celerity_control_fails_only_the_new_criterion(probe, tmp_path):
    outcome = run_probe(
        registry.find_model("reference_fixed_celerity"),
        probe,
        gate_seeds(PROBE_ID, 1),
        workdir=tmp_path,
    )
    assert outcome.verdict == FAIL
    assert [item.name for item in outcome.failing] == ["wave_celerity_bounds"]


def test_saint_venant_passes_one_gate_seed_end_to_end(probe, tmp_path):
    outcome = run_probe(
        registry.find_model("reference_saint_venant"),
        probe,
        gate_seeds(PROBE_ID, 1),
        workdir=tmp_path,
    )
    assert outcome.verdict == PASS, outcome.failing
