"""Unit and integration tests for FieldInfluenceAnalysis (Etap A3).

Strategy
--------
Full VQE runs are slow (seconds each) and non-deterministic.  The tests here
use two approaches:

1. **Unit tests** - test each internal helper in isolation using mocking or
   small, purely-classical computations (no VQE).
2. **Smoke / integration test** - one short end-to-end run with a tiny chain
   and very few VQE iterations to verify the full pipeline glues together.

The smoke test is marked ``@pytest.mark.slow`` so it can be excluded from fast
CI runs::

    pytest tests/ -m "not slow"   # skip slow
    pytest tests/ -m slow         # only slow
"""

from __future__ import annotations

import dataclasses
import math
import os
from typing import Any
from unittest.mock import patch

import pytest

# Force a non-interactive matplotlib backend before any matplotlib import so
# that plot tests work in headless environments (no Tk / display server).
os.environ.setdefault("MPLBACKEND", "Agg")

# NOTE: Use non-src-prefixed imports (matching pytest pythonpath = ["src"]).
# Mixing src.* and bare imports produces two separate class objects in Python's
# module cache, which breaks isinstance() checks.
from analysis.field_influence_analysis import (
    FieldInfluenceAnalysis,
    ScenarioResult,
)
from enums import InteractionType
from particle.external_field import ExternalField

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_SHORT_CHAIN = "HPPHH"  # length-5 HP chain - minimum for backbone contacts
_SIDE_CHAIN = "_____"


def _make_analysis(**kwargs) -> FieldInfluenceAnalysis:
    """Construct a FieldInfluenceAnalysis with sensible test defaults."""
    defaults = {
        "main_chain": _SHORT_CHAIN,
        "interaction_type": InteractionType.HP,
        "uniform_lambdas": [0.5, 1.0],
        "vqe_max_iter": 5,
    }
    defaults.update(kwargs)
    return FieldInfluenceAnalysis(**defaults)


def _fake_scenario_result(label: str, energy: float = -1.0) -> ScenarioResult:
    """Build a minimal ScenarioResult for plotting tests."""
    return ScenarioResult(
        label=label,
        field=None,
        minimum_energy=energy,
        best_bitstring="0101",
        state_probabilities={"0101": 0.7, "1010": 0.2, "0011": 0.1},
        vqe_iterations=[1, 2, 3],
        vqe_energies=[-0.5, -0.8, energy],
    )


# ===========================================================================
# __init__ - construction
# ===========================================================================


class TestInit:
    def test_default_uniform_lambdas(self):
        a = FieldInfluenceAnalysis(main_chain=_SHORT_CHAIN)
        assert a.uniform_lambdas == [0.1, 0.5, 1.0, 2.0]

    def test_custom_uniform_lambdas(self):
        a = _make_analysis(uniform_lambdas=[0.5, 2.0])
        assert a.uniform_lambdas == [0.5, 2.0]

    def test_default_side_chain_all_placeholders(self):
        a = _make_analysis()
        assert a.side_chain == "_" * len(_SHORT_CHAIN)

    def test_custom_side_chain(self):
        a = FieldInfluenceAnalysis(main_chain=_SHORT_CHAIN, side_chain=_SIDE_CHAIN)
        assert a.side_chain == _SIDE_CHAIN

    def test_interaction_type_stored(self):
        a = _make_analysis(interaction_type=InteractionType.HP)
        assert a.interaction_type == InteractionType.HP

    def test_results_empty_before_run(self):
        a = _make_analysis()
        assert a.results == []

    def test_protein_chain_length_matches(self):
        a = _make_analysis()
        assert len(a._protein.main_chain) == len(_SHORT_CHAIN)

    def test_mj_interaction_type_accepted(self):
        # MJ requires valid MJ symbols; use a short sequence known to be valid.
        a = FieldInfluenceAnalysis(
            main_chain="AAAAA",
            interaction_type=InteractionType.MJ,
            uniform_lambdas=[0.5],
            vqe_max_iter=2,
        )
        assert a.interaction_type == InteractionType.MJ


# ===========================================================================
# _build_non_uniform_field
# ===========================================================================


class TestBuildNonUniformField:
    def test_returns_external_field(self):
        a = _make_analysis()
        field = a._build_non_uniform_field()
        assert isinstance(field, ExternalField)

    def test_has_correct_number_of_nodes(self):
        a = _make_analysis()
        field = a._build_non_uniform_field()
        # One node per bead in the main chain
        assert len(field.nodes()) == len(_SHORT_CHAIN)

    def test_central_bead_has_most_negative_energy(self):
        """The mid-chain bead should have the most negative (strongest) field energy."""
        a = _make_analysis()
        field = a._build_non_uniform_field()
        chain_len = len(_SHORT_CHAIN)
        mid = (chain_len - 1) // 2
        energies = {i: field.get_energy((i,)) for i in range(chain_len)}
        centre_e = energies[mid]
        for i, e in energies.items():
            assert e >= centre_e - 1e-9, (
                f"Bead {i} (e={e}) is more negative than centre bead {mid} (e={centre_e})"
            )

    def test_all_energies_finite(self):
        a = _make_analysis()
        field = a._build_non_uniform_field()
        for i in range(len(_SHORT_CHAIN)):
            e = field.get_energy((i,))
            assert math.isfinite(e), f"Bead {i} has non-finite field energy: {e}"

    def test_terminal_beads_weaker_than_centre(self):
        """Termini (index 0 and N-1) should be weaker than the centre bead."""
        a = _make_analysis()
        field = a._build_non_uniform_field()
        n = len(_SHORT_CHAIN)
        mid = (n - 1) // 2
        e_mid = field.get_energy((mid,))
        e_term_left = field.get_energy((0,))
        e_term_right = field.get_energy((n - 1,))
        assert e_term_left > e_mid, (
            "Left terminus should be weaker (less negative) than centre"
        )
        assert e_term_right > e_mid, (
            "Right terminus should be weaker (less negative) than centre"
        )


# ===========================================================================
# summary()  # noqa: ERA001
# ===========================================================================


class TestSummary:
    def test_summary_raises_before_run(self):
        a = _make_analysis()
        with pytest.raises(RuntimeError, match=r"run\(\)"):
            a.summary()

    def test_summary_contains_all_labels(self):
        a = _make_analysis()
        a.results = [
            _fake_scenario_result("baseline (no field)", -1.0),
            _fake_scenario_result("λ=0.5 (uniform)", -1.3),
            _fake_scenario_result("non-uniform (centre boost)", -1.1),
        ]
        s = a.summary()
        assert "baseline (no field)" in s
        assert "λ=0.5 (uniform)" in s
        assert "non-uniform (centre boost)" in s

    def test_summary_contains_energies(self):
        a = _make_analysis()
        a.results = [_fake_scenario_result("test", -2.345678)]
        s = a.summary()
        assert "-2.345678" in s


# ===========================================================================
# plot() - guard against calling before run()
# ===========================================================================


class TestPlotGuard:
    def test_plot_raises_before_run(self):
        a = _make_analysis()
        with pytest.raises(RuntimeError, match=r"run\(\)"):
            a.plot()


# ===========================================================================
# plot() - output saved to disk when output_dir is provided
# ===========================================================================


class TestPlotOutput:
    """Tests that plot() actually writes PNG files to disk."""

    def _populate_results(self, analysis: FieldInfluenceAnalysis) -> None:
        analysis.results = [
            _fake_scenario_result("baseline (no field)", -1.0),
            _fake_scenario_result("λ=0.5 (uniform)", -1.2),
            _fake_scenario_result("λ=1.0 (uniform)", -1.4),
            _fake_scenario_result("non-uniform (centre boost)", -1.1),
        ]

    def test_energy_vs_lambda_png_created(self, tmp_path):
        a = _make_analysis()
        self._populate_results(a)
        a.plot(output_dir=tmp_path)
        assert (tmp_path / "energy_vs_lambda.png").exists()

    def test_energy_bar_png_created(self, tmp_path):
        a = _make_analysis()
        self._populate_results(a)
        a.plot(output_dir=tmp_path)
        assert (tmp_path / "energy_comparison_bar.png").exists()

    def test_probability_distributions_png_created(self, tmp_path):
        a = _make_analysis()
        self._populate_results(a)
        a.plot(output_dir=tmp_path)
        assert (tmp_path / "probability_distributions.png").exists()

    def test_output_dir_created_if_missing(self, tmp_path):
        a = _make_analysis()
        self._populate_results(a)
        new_dir = tmp_path / "new" / "nested" / "dir"
        assert not new_dir.exists()
        a.plot(output_dir=new_dir)
        assert new_dir.exists()

    def test_only_uniform_results_still_plots(self, tmp_path):
        """If there are no non-uniform results the plot should still work."""
        a = _make_analysis()
        a.results = [
            _fake_scenario_result("baseline (no field)", -1.0),
            _fake_scenario_result("λ=0.5 (uniform)", -1.2),
        ]
        a.plot(output_dir=tmp_path)
        assert (tmp_path / "energy_vs_lambda.png").exists()


# ===========================================================================
# ScenarioResult dataclass
# ===========================================================================


class TestScenarioResult:
    def test_is_dataclass(self):
        assert dataclasses.is_dataclass(ScenarioResult)

    def test_fields_accessible(self):
        r = _fake_scenario_result("test", -0.5)
        assert r.label == "test"
        assert r.minimum_energy == pytest.approx(-0.5)
        assert r.best_bitstring == "0101"
        assert isinstance(r.state_probabilities, dict)
        assert isinstance(r.vqe_iterations, list)
        assert isinstance(r.vqe_energies, list)

    def test_field_can_be_none(self):
        r = _fake_scenario_result("baseline")
        assert r.field is None

    def test_field_can_be_external_field(self):
        field = ExternalField.uniform(strength=-1.0)
        r = ScenarioResult(
            label="test",
            field=field,
            minimum_energy=-1.0,
            best_bitstring="01",
            state_probabilities={"01": 1.0},
            vqe_iterations=[1],
            vqe_energies=[-1.0],
        )
        assert r.field is field


# ===========================================================================
# run() - mocked VQE to verify orchestration logic
# ===========================================================================


class TestRunOrchestration:
    """Verify that run() calls _run_scenario the right number of times
    and with the right arguments, without actually executing VQE."""

    def _mock_scenario_result(self, label: str, field: Any) -> ScenarioResult:
        return _fake_scenario_result(label)

    def test_number_of_scenarios_called(self):
        lambdas = [0.1, 0.5, 1.0]
        a = _make_analysis(uniform_lambdas=lambdas)

        with patch.object(a, "_run_scenario", side_effect=self._mock_scenario_result):
            a.run()

        # 1 baseline + len(lambdas) uniform + 1 non-uniform
        assert len(a.results) == 1 + len(lambdas) + 1

    def test_baseline_is_first(self):
        a = _make_analysis(uniform_lambdas=[1.0])
        with patch.object(a, "_run_scenario", side_effect=self._mock_scenario_result):
            a.run()
        assert a.results[0].label == "baseline (no field)"

    def test_non_uniform_is_last(self):
        a = _make_analysis(uniform_lambdas=[1.0])
        with patch.object(a, "_run_scenario", side_effect=self._mock_scenario_result):
            a.run()
        assert "non-uniform" in a.results[-1].label

    def test_uniform_labels_contain_lambda_values(self):
        lambdas = [0.5, 2.0]
        a = _make_analysis(uniform_lambdas=lambdas)
        with patch.object(a, "_run_scenario", side_effect=self._mock_scenario_result):
            a.run()
        # results[1] and results[2] should be the uniform scenarios
        for r, lam in zip(a.results[1:-1], lambdas, strict=False):
            assert str(lam) in r.label

    def test_run_clears_previous_results(self):
        a = _make_analysis(uniform_lambdas=[1.0])
        a.results = [_fake_scenario_result("old")]
        with patch.object(a, "_run_scenario", side_effect=self._mock_scenario_result):
            a.run()
        assert all(r.label != "old" for r in a.results)

    def test_baseline_field_arg_is_none(self):
        """_run_scenario for the baseline must receive field=None."""
        a = _make_analysis(uniform_lambdas=[])
        captured_fields: list[Any] = []

        def _capture(label: str, field: Any) -> ScenarioResult:
            captured_fields.append(field)
            return _fake_scenario_result(label)

        with patch.object(a, "_run_scenario", side_effect=_capture):
            a.run()

        # First call is baseline, last is non-uniform; no uniform ones
        assert captured_fields[0] is None

    def test_uniform_scenarios_receive_external_field(self):
        lambdas = [0.5, 1.0]
        a = _make_analysis(uniform_lambdas=lambdas)
        captured_fields: list[Any] = []

        def _capture(label: str, field: Any) -> ScenarioResult:
            captured_fields.append(field)
            return _fake_scenario_result(label)

        with patch.object(a, "_run_scenario", side_effect=_capture):
            a.run()

        # indices 1 and 2 are the uniform scenarios
        for field in captured_fields[1:-1]:
            assert isinstance(field, ExternalField), (
                f"Expected ExternalField for uniform scenario, got {type(field)}"
            )

    def test_non_uniform_scenario_receives_external_field(self):
        a = _make_analysis(uniform_lambdas=[])
        captured_fields: list[Any] = []

        def _capture(label: str, field: Any) -> ScenarioResult:
            captured_fields.append(field)
            return _fake_scenario_result(label)

        with patch.object(a, "_run_scenario", side_effect=_capture):
            a.run()

        # Last call is non-uniform
        assert isinstance(captured_fields[-1], ExternalField)


# ===========================================================================
# Smoke / integration test (slow - real VQE, few iterations)
# ===========================================================================


@pytest.mark.slow
class TestEndToEnd:
    """Full pipeline smoke test: construct → run (5 VQE iterations) → plot."""

    def test_run_and_plot(self, tmp_path):
        analysis = FieldInfluenceAnalysis(
            main_chain="HPPHH",
            interaction_type=InteractionType.HP,
            uniform_lambdas=[0.5, 1.0],
            vqe_max_iter=5,
        )
        analysis.run()

        # Basic sanity on results structure
        assert len(analysis.results) == 1 + 2 + 1  # baseline + 2 uniform + non-uniform
        for r in analysis.results:
            assert isinstance(r.minimum_energy, float)
            assert math.isfinite(r.minimum_energy)
            assert isinstance(r.best_bitstring, str)

        # Plots must be produced without error
        analysis.plot(output_dir=tmp_path)
        assert (tmp_path / "energy_vs_lambda.png").exists()
        assert (tmp_path / "energy_comparison_bar.png").exists()
        assert (tmp_path / "probability_distributions.png").exists()

        # Summary must not raise
        s = analysis.summary()
        assert "baseline" in s
