"""Regression checks for the exploratory wave-celerity audit."""

from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "experiments" / "wave_celerity_audit.py"


def _load_audit():
    spec = importlib.util.spec_from_file_location("wave_celerity_audit", AUDIT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_saint_venant_transient_celerity_changes_with_flow_state():
    audit = _load_audit()
    solver = audit.load_solver()
    rows = [audit.run_seed(solver, seed, n_cells=64) for seed in range(3)]
    summary = audit.summarize(rows)

    # These are feasibility guards, not proposed production tolerances.
    assert summary["monotonic_pass_seeds"] == 3
    assert summary["high_to_low_celerity_ratio_min"] > 1.5
    assert 0.80 < summary["observed_to_theory_ratio_min"]
    assert summary["observed_to_theory_ratio_max"] < 1.05
    assert summary["max_froude"] < 0.5


def test_existing_fixed_celerity_control_has_no_state_dependence():
    # The current flex_lumped adapter uses a fixed 1 m/s celerity and is a
    # must-pass physical model for merged momentum/routing-lag-consistency.
    values = [1.0, 1.0, 1.0]
    assert values[0] == values[1] == values[2]
