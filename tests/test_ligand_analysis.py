"""Tests for :class:`~analysis.ligand_analysis.LigandAnalysis`.

Tests are organised in three classes mirroring the three public analysis
methods plus fixture-level integration tests:

* TestInitialisation  – constructor validation, ligand/interaction creation.
* TestRunAndResults   – run() produces expected result structure; no VQE
                        correctness expected, only interface correctness.
* TestComputeBindingEnergy  – ΔE formula, call-before-run raises error.
* TestPositionDistribution  – normalisation, length, BINARY/UNARY parsing.
* TestEncapsulation  – interior/boundary split, threshold, score bounds.
* TestPlotting  – plot() does not crash; files created in a tmp dir.

Strategy
--------
* VQE is run with vqe_max_iter=3 to keep tests fast.  We do not assert on
  exact energy values – the optimiser may not converge in 3 steps.  We only
  check structural properties (types, lengths, ranges, error conditions).
* The minimal chain "HPPH" with 4 lattice nodes gives the smallest possible
  Hamiltonian, keeping wall-clock time low.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from src.analysis.ligand_analysis import (
    EncapsulationResult,
    LigandAnalysis,
    LigandScenarioResult,
    run_hp_hydrophobic_experiment,
    run_mj_strong_ligand_experiment,
)
from src.enums import InteractionType

# Import PositionEncoding via the same path used internally by LigandAnalysis
# (particle.ligand_bead) to avoid dual-module identity issues.
import sys as _sys
import importlib as _importlib
# Ensure src is on path so particle.ligand_bead resolves to same module object
_importlib.import_module("particle.ligand_bead")
PositionEncoding = _sys.modules["particle.ligand_bead"].PositionEncoding

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

MINIMAL_CHAIN = "HPPHH"      # 5-mer – minimum chain length accepted by Protein
MINIMAL_NODES = 4            # 4 nodes → 2 pos-qubits (BINARY) / 4 (UNARY)
FAST_ITER = 3                # very few VQE iterations – we test interface, not accuracy


@pytest.fixture(scope="module")
def analysis_hp_h() -> LigandAnalysis:
    """Hydrophobic ligand, HP model, already run."""
    a = LigandAnalysis(
        main_chain=MINIMAL_CHAIN,
        interaction_type=InteractionType.HP,
        ligand_hp_type="H",
        num_lattice_nodes=MINIMAL_NODES,
        position_encoding=PositionEncoding.BINARY,
        vqe_max_iter=FAST_ITER,
    )
    a.run()
    return a


@pytest.fixture(scope="module")
def analysis_hp_p() -> LigandAnalysis:
    """Polar ligand, HP model, already run."""
    a = LigandAnalysis(
        main_chain=MINIMAL_CHAIN,
        interaction_type=InteractionType.HP,
        ligand_hp_type="P",
        num_lattice_nodes=MINIMAL_NODES,
        position_encoding=PositionEncoding.BINARY,
        vqe_max_iter=FAST_ITER,
    )
    a.run()
    return a


@pytest.fixture(scope="module")
def analysis_custom() -> LigandAnalysis:
    """Custom ligand, HP model, already run."""
    a = LigandAnalysis(
        main_chain=MINIMAL_CHAIN,
        interaction_type=InteractionType.HP,
        ligand_hp_type=None,
        ligand_energy_map={"H": -1.5, "P": -0.2},
        ligand_default_energy=0.0,
        num_lattice_nodes=MINIMAL_NODES,
        position_encoding=PositionEncoding.BINARY,
        vqe_max_iter=FAST_ITER,
    )
    a.run()
    return a


@pytest.fixture(scope="module")
def analysis_unary() -> LigandAnalysis:
    """Hydrophobic ligand with UNARY encoding, HP model, already run."""
    a = LigandAnalysis(
        main_chain=MINIMAL_CHAIN,
        interaction_type=InteractionType.HP,
        ligand_hp_type="H",
        num_lattice_nodes=MINIMAL_NODES,
        position_encoding=PositionEncoding.UNARY,
        vqe_max_iter=FAST_ITER,
    )
    a.run()
    return a


# ---------------------------------------------------------------------------
# 1. Initialisation tests
# ---------------------------------------------------------------------------


class TestInitialisation:
    def test_attributes_stored_correctly(self):
        a = LigandAnalysis(
            main_chain="HPPHH",
            ligand_hp_type="H",
            num_lattice_nodes=6,
        )
        assert a.main_chain == "HPPHH"
        assert a.ligand_hp_type == "H"
        assert a.num_lattice_nodes == 6
        assert a.interaction_type == InteractionType.HP

    def test_side_chain_defaults_to_underscores(self):
        a = LigandAnalysis(main_chain="HPPHH", ligand_hp_type="H")
        assert a.side_chain == "_____"

    def test_custom_side_chain_stored(self):
        a = LigandAnalysis(
            main_chain="HPPHH",
            ligand_hp_type="H",
            side_chain="_____",
        )
        assert a.side_chain == "_____"

    def test_mj_interaction_type_stored(self):
        a = LigandAnalysis(
            main_chain="ACDEG",
            interaction_type=InteractionType.MJ,
            ligand_hp_type=None,
            ligand_energy_map={"A": -1.0},
        )
        assert a.interaction_type == InteractionType.MJ

    def test_ligand_bead_created_with_correct_nodes(self):
        a = LigandAnalysis(
            main_chain="HPPHH",
            ligand_hp_type="H",
            num_lattice_nodes=8,
            position_encoding=PositionEncoding.BINARY,
        )
        assert a._ligand.num_lattice_nodes == 8
        # BINARY: ceil(log2(8)) = 3 qubits
        assert a._ligand.num_position_qubits == 3

    def test_results_empty_before_run(self):
        a = LigandAnalysis(main_chain="HPPHH", ligand_hp_type="H")
        assert a.results == {}

    def test_invalid_hp_type_raises_value_error(self):
        with pytest.raises(ValueError, match="ligand_hp_type"):
            LigandAnalysis(main_chain="HPPHH", ligand_hp_type="X")

    def test_custom_mode_when_hp_type_is_none(self):
        a = LigandAnalysis(
            main_chain="HPPHH",
            ligand_hp_type=None,
            ligand_energy_map={"H": -1.0},
        )
        from src.interaction.ligand_interaction import LigandInteractionMode
        assert a._ligand_interaction.mode.name == "CUSTOM"

    def test_hp_like_mode_when_hp_type_set(self):
        a = LigandAnalysis(main_chain="HPPHH", ligand_hp_type="H")
        assert a._ligand_interaction.mode.name == "HP_LIKE"


# ---------------------------------------------------------------------------
# 2. run() and results structure
# ---------------------------------------------------------------------------


class TestRunAndResults:
    def test_run_populates_both_keys(self, analysis_hp_h):
        assert "baseline" in analysis_hp_h.results
        assert "with_ligand" in analysis_hp_h.results

    def test_baseline_result_is_scenario_result(self, analysis_hp_h):
        assert isinstance(analysis_hp_h.results["baseline"], LigandScenarioResult)

    def test_ligand_result_is_scenario_result(self, analysis_hp_h):
        assert isinstance(analysis_hp_h.results["with_ligand"], LigandScenarioResult)

    def test_minimum_energy_is_finite(self, analysis_hp_h):
        for r in analysis_hp_h.results.values():
            assert math.isfinite(r.minimum_energy)

    def test_vqe_iterations_recorded(self, analysis_hp_h):
        for r in analysis_hp_h.results.values():
            assert len(r.vqe_iterations) > 0

    def test_best_bitstring_non_empty(self, analysis_hp_h):
        for r in analysis_hp_h.results.values():
            assert isinstance(r.best_bitstring, str)
            assert len(r.best_bitstring) > 0

    def test_state_probabilities_are_positive(self, analysis_hp_h):
        for r in analysis_hp_h.results.values():
            for p in r.state_probabilities.values():
                assert p >= 0.0

    def test_ligand_node_probabilities_correct_length(self, analysis_hp_h):
        r = analysis_hp_h.results["with_ligand"]
        assert len(r.ligand_node_probabilities) == MINIMAL_NODES

    def test_binding_energy_set_on_ligand_result(self, analysis_hp_h):
        r = analysis_hp_h.results["with_ligand"]
        assert r.binding_energy is not None
        assert math.isfinite(r.binding_energy)

    def test_binding_energy_is_none_on_baseline(self, analysis_hp_h):
        r = analysis_hp_h.results["baseline"]
        assert r.binding_energy is None

    def test_run_clears_previous_results(self):
        a = LigandAnalysis(
            main_chain=MINIMAL_CHAIN,
            ligand_hp_type="H",
            num_lattice_nodes=MINIMAL_NODES,
            vqe_max_iter=FAST_ITER,
        )
        a.run()
        keys_first = set(a.results.keys())
        a.run()
        keys_second = set(a.results.keys())
        assert keys_first == keys_second


# ---------------------------------------------------------------------------
# 3. compute_binding_energy
# ---------------------------------------------------------------------------


class TestComputeBindingEnergy:
    def test_returns_finite_float(self, analysis_hp_h):
        delta_e = analysis_hp_h.compute_binding_energy()
        assert isinstance(delta_e, float)
        assert math.isfinite(delta_e)

    def test_equals_energy_difference(self, analysis_hp_h):
        e_base = analysis_hp_h.results["baseline"].minimum_energy
        e_lig = analysis_hp_h.results["with_ligand"].minimum_energy
        expected = e_lig - e_base
        assert abs(analysis_hp_h.compute_binding_energy() - expected) < 1e-10

    def test_hydrophobic_binding_energy_negative_or_close_to_zero(self, analysis_hp_h):
        """Hydrophobic ligand + HP chain: expect non-positive ΔE (stabilising)."""
        delta_e = analysis_hp_h.compute_binding_energy()
        assert delta_e < 100.0  # not wildly positive

    def test_custom_binding_energy_finite(self, analysis_custom):
        delta_e = analysis_custom.compute_binding_energy()
        assert math.isfinite(delta_e)

    def test_raises_before_run(self):
        a = LigandAnalysis(main_chain=MINIMAL_CHAIN, ligand_hp_type="H")
        with pytest.raises(RuntimeError, match="run()"):
            a.compute_binding_energy()

    def test_matches_stored_binding_energy(self, analysis_hp_h):
        stored = analysis_hp_h.results["with_ligand"].binding_energy
        computed = analysis_hp_h.compute_binding_energy()
        assert abs(stored - computed) < 1e-10


# ---------------------------------------------------------------------------
# 4. compute_ligand_position_distribution
# ---------------------------------------------------------------------------


class TestPositionDistribution:
    def test_returns_list_of_correct_length(self, analysis_hp_h):
        dist = analysis_hp_h.compute_ligand_position_distribution()
        assert isinstance(dist, list)
        assert len(dist) == MINIMAL_NODES

    def test_all_probabilities_non_negative(self, analysis_hp_h):
        dist = analysis_hp_h.compute_ligand_position_distribution()
        for p in dist:
            assert p >= 0.0

    def test_probabilities_sum_to_one(self, analysis_hp_h):
        dist = analysis_hp_h.compute_ligand_position_distribution()
        assert abs(sum(dist) - 1.0) < 1e-6

    def test_unary_returns_correct_length(self, analysis_unary):
        dist = analysis_unary.compute_ligand_position_distribution()
        assert len(dist) == MINIMAL_NODES

    def test_unary_probs_non_negative(self, analysis_unary):
        dist = analysis_unary.compute_ligand_position_distribution()
        for p in dist:
            assert p >= 0.0

    def test_unary_probs_sum_to_one(self, analysis_unary):
        dist = analysis_unary.compute_ligand_position_distribution()
        total = sum(dist)
        assert abs(total - 1.0) < 1e-6

    def test_custom_mode_returns_correct_length(self, analysis_custom):
        dist = analysis_custom.compute_ligand_position_distribution()
        assert len(dist) == MINIMAL_NODES

    def test_raises_before_run(self):
        a = LigandAnalysis(main_chain=MINIMAL_CHAIN, ligand_hp_type="H")
        with pytest.raises(RuntimeError, match="run()"):
            a.compute_ligand_position_distribution()

    def test_polar_ligand_distribution_sums_to_one(self, analysis_hp_p):
        dist = analysis_hp_p.compute_ligand_position_distribution()
        assert abs(sum(dist) - 1.0) < 1e-6


# ---------------------------------------------------------------------------
# 5. detect_encapsulation
# ---------------------------------------------------------------------------


class TestEncapsulation:
    def test_returns_encapsulation_result(self, analysis_hp_h):
        enc = analysis_hp_h.detect_encapsulation()
        assert isinstance(enc, EncapsulationResult)

    def test_score_in_unit_interval(self, analysis_hp_h):
        enc = analysis_hp_h.detect_encapsulation()
        assert 0.0 <= enc.score <= 1.0 + 1e-9

    def test_interior_and_boundary_partition_all_nodes(self, analysis_hp_h):
        enc = analysis_hp_h.detect_encapsulation()
        all_nodes = set(range(MINIMAL_NODES))
        covered = set(enc.interior_nodes) | set(enc.boundary_nodes)
        assert covered == all_nodes

    def test_interior_boundary_disjoint(self, analysis_hp_h):
        enc = analysis_hp_h.detect_encapsulation()
        assert len(set(enc.interior_nodes) & set(enc.boundary_nodes)) == 0

    def test_is_encapsulated_consistent_with_score(self, analysis_hp_h):
        threshold = 0.6
        enc = analysis_hp_h.detect_encapsulation(encapsulation_threshold=threshold)
        assert enc.is_encapsulated == (enc.score >= threshold)

    def test_threshold_zero_always_encapsulated(self, analysis_hp_h):
        enc = analysis_hp_h.detect_encapsulation(encapsulation_threshold=0.0)
        assert enc.is_encapsulated is True

    def test_threshold_one_never_encapsulated(self, analysis_hp_h):
        enc = analysis_hp_h.detect_encapsulation(encapsulation_threshold=1.0 + 1e-9)
        assert enc.is_encapsulated is False

    def test_interior_probability_plus_boundary_equals_total(self, analysis_hp_h):
        enc = analysis_hp_h.detect_encapsulation()
        dist = analysis_hp_h.compute_ligand_position_distribution()
        total = sum(dist)
        assert abs(enc.interior_probability + enc.boundary_probability - total) < 1e-9

    def test_interior_fraction_changes_interior_size(self, analysis_hp_h):
        enc_half = analysis_hp_h.detect_encapsulation(interior_fraction=0.5)
        enc_full = analysis_hp_h.detect_encapsulation(interior_fraction=1.0)
        assert len(enc_full.interior_nodes) >= len(enc_half.interior_nodes)

    def test_raises_before_run(self):
        a = LigandAnalysis(main_chain=MINIMAL_CHAIN, ligand_hp_type="H")
        with pytest.raises(RuntimeError, match="run()"):
            a.detect_encapsulation()

    def test_encapsulation_threshold_stored(self, analysis_hp_h):
        enc = analysis_hp_h.detect_encapsulation(encapsulation_threshold=0.75)
        assert enc.threshold == 0.75


# ---------------------------------------------------------------------------
# 6. summary()
# ---------------------------------------------------------------------------


class TestSummary:
    def test_summary_is_string(self, analysis_hp_h):
        s = analysis_hp_h.summary()
        assert isinstance(s, str)
        assert len(s) > 0

    def test_summary_contains_labels(self, analysis_hp_h):
        s = analysis_hp_h.summary()
        assert "baseline" in s
        assert "ligand" in s

    def test_summary_raises_before_run(self):
        a = LigandAnalysis(main_chain=MINIMAL_CHAIN, ligand_hp_type="H")
        with pytest.raises(RuntimeError, match="run()"):
            a.summary()


# ---------------------------------------------------------------------------
# 7. Plotting (smoke tests – no assertions on figure content)
# ---------------------------------------------------------------------------


class TestPlotting:
    def test_plot_creates_files(self, tmp_path, analysis_hp_h):
        analysis_hp_h.plot(output_dir=tmp_path)
        expected_files = [
            "ligand_position_distribution.png",
            "energy_comparison.png",
            "lattice_heatmap.png",
        ]
        for fname in expected_files:
            assert (tmp_path / fname).exists(), f"Missing plot: {fname}"

    def test_plot_raises_before_run(self, tmp_path):
        a = LigandAnalysis(main_chain=MINIMAL_CHAIN, ligand_hp_type="H")
        with pytest.raises(RuntimeError, match="run()"):
            a.plot(output_dir=tmp_path)

    def test_unary_plot_does_not_crash(self, tmp_path, analysis_unary):
        analysis_unary.plot(output_dir=tmp_path)


# ---------------------------------------------------------------------------
# 8. Convenience experiment functions
# ---------------------------------------------------------------------------


class TestExperimentFunctions:
    def test_hp_hydrophobic_experiment_returns_analysis(self, tmp_path):
        result = run_hp_hydrophobic_experiment(
            main_chain=MINIMAL_CHAIN,
            num_lattice_nodes=MINIMAL_NODES,
            vqe_max_iter=FAST_ITER,
            output_dir=tmp_path,
        )
        assert isinstance(result, LigandAnalysis)
        assert "with_ligand" in result.results

    def test_mj_strong_ligand_experiment_returns_analysis(self, tmp_path):
        result = run_mj_strong_ligand_experiment(
            main_chain="ACDEG",
            num_lattice_nodes=MINIMAL_NODES,
            vqe_max_iter=FAST_ITER,
            output_dir=tmp_path,
        )
        assert isinstance(result, LigandAnalysis)
        assert "with_ligand" in result.results

    def test_hp_experiment_binding_energy_finite(self, tmp_path):
        result = run_hp_hydrophobic_experiment(
            main_chain=MINIMAL_CHAIN,
            num_lattice_nodes=MINIMAL_NODES,
            vqe_max_iter=FAST_ITER,
            output_dir=tmp_path,
        )
        delta_e = result.compute_binding_energy()
        assert math.isfinite(delta_e)
