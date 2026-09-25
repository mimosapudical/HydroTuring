"""Focused tests for transient wave-celerity bounds."""

from __future__ import annotations

import importlib.util

import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria import get, is_paired
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


def _load_adapter(model_name: str):
    path = registry.find_model(model_name).path / "ht_adapter.py"
    spec = importlib.util.spec_from_file_location(f"{model_name}_module", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_criterion_is_paired_and_is_the_probe_headline(probe):
    assert is_paired("wave_celerity_bounds")
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
def test_short_long_pairs_are_identical_except_reach_length(probe, seed):
    for state in STATES:
        short = build_case(probe, seed, f"{state}_short")
        long = build_case(probe, seed, f"{state}_long")
        pd.testing.assert_frame_equal(short.forcing, long.forcing, check_exact=True)
        assert set(short.static) == set(long.static)
        for key in short.static:
            if key != "reach_length_m":
                assert short.static[key] == long.static[key]
        assert short.static["reach_length_m"] == pytest.approx(10_000.0)
        assert long.static["reach_length_m"] == pytest.approx(50_000.0)


@pytest.mark.parametrize("seed", gate_seeds(PROBE_ID, 3))
def test_operating_states_raise_only_the_base_hydraulic_forcing(probe, seed):
    cases = [build_case(probe, seed, f"{state}_short") for state in STATES]
    effective = [
        (case.forcing["pr"] - case.forcing["pet"]).iloc[case.spinup_steps + 1]
        for case in cases
    ]
    assert effective[1] == pytest.approx(2.0 * effective[0])
    assert effective[2] == pytest.approx(4.0 * effective[0])
    for case in cases:
        scored = case.after_spinup(case.forcing)
        assert int((scored["_pulse"] > 0.0).sum()) == 2
        assert float(scored["_pulse"].max()) == pytest.approx(
            0.05 * float((scored["pr"] - scored["pet"]).min())
        )


def test_only_process_aware_references_are_compatible(probe):
    case = build_case(probe, gate_seeds(PROBE_ID, 1)[0], "low_short")
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
    discharge = 20.0
    exact = module._rectangular_celerity(discharge, depth, width)

    # Differentiate Q(A) numerically while holding n and S fixed. The common
    # Manning prefactor cancels after calibrating it at the chosen state.
    area = width * depth
    perimeter = width + 2.0 * depth
    radius = area / perimeter
    factor = discharge / (area * radius ** (2.0 / 3.0))

    def rating(a):
        d = a / width
        r = a / (width + 2.0 * d)
        return factor * a * r ** (2.0 / 3.0)

    eps = 1.0e-4 * area
    finite_difference = (rating(area + eps) - rating(area - eps)) / (2.0 * eps)
    assert exact == pytest.approx(finite_difference, rel=1.0e-8)


def test_fixed_celerity_control_fails_the_registered_probe(probe, tmp_path):
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
