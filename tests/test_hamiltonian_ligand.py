"""Unit tests for the LigandBead ↔ HamiltonianBuilder integration (Etap B3).

Tests verify:
- sum_hamiltonians() without ligand is backward compatible (same qubit count,
  same operator as without ligand argument).
- sum_hamiltonians(ligand, ligand_interaction) correctly increases the qubit
  count by ligand.num_position_qubits.
- _build_ligand_contact_term() returns a valid SparsePauliOp.
- _build_ligand_contact_term() with no ligand returns a zero operator of the
  protein-only qubit size.
- _build_ligand_contact_term() with a ligand returns an operator acting on
  n_turn + n_position qubits.
- The identity coefficient of H_ligand equals sum_i E_i (Variant A property).
- Passing only ligand or only ligand_interaction raises ValueError.
- Both HP-like and Custom ligand interactions are tested.
- Qubit count formula: n_total = (N-1)*QUBITS_PER_TURN + ligand.num_position_qubits.
"""

from __future__ import annotations

import pytest
from qiskit.quantum_info import SparsePauliOp

from src.builder.hamiltonian_builder import HamiltonianBuilder
from src.constants import QUBITS_PER_TURN
from src.contact.contact_map import ContactMap
from src.distance.distance_map import DistanceMap
from src.interaction.hp_interaction import HPInteraction
from src.interaction.ligand_interaction import LigandInteraction
from src.particle.ligand_bead import LigandBead, PositionEncoding
from src.protein.protein import Protein
from src.utils.qubit_utils import pad_to_n_qubits

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

MAIN_CHAIN = "HPPHH"  # length 5 - minimum to allow backbone contacts
SIDE_CHAIN = "_____"
NUM_LATTICE_NODES = 4  # keeps ligand qubit count small (2 for BINARY, 4 for UNARY)


@pytest.fixture(scope="module")
def hp_interaction() -> HPInteraction:
    return HPInteraction()


@pytest.fixture(scope="module")
def protein(hp_interaction: HPInteraction) -> Protein:
    return Protein(
        main_protein_sequence=MAIN_CHAIN,
        side_protein_sequence=SIDE_CHAIN,
        valid_symbols=hp_interaction.valid_symbols,
    )


@pytest.fixture(scope="module")
def contact_map(protein: Protein) -> ContactMap:
    return ContactMap(protein=protein)


@pytest.fixture(scope="module")
def distance_map(protein: Protein) -> DistanceMap:
    return DistanceMap(protein=protein)


@pytest.fixture(scope="module")
def ligand_binary() -> LigandBead:
    """Ligand with BINARY encoding - 2 position qubits for 4 nodes."""
    return LigandBead(
        symbol="L",
        index=0,
        num_lattice_nodes=NUM_LATTICE_NODES,
        encoding=PositionEncoding.BINARY,
    )


@pytest.fixture(scope="module")
def ligand_unary() -> LigandBead:
    """Ligand with UNARY encoding - 4 position qubits for 4 nodes."""
    return LigandBead(
        symbol="L",
        index=0,
        num_lattice_nodes=NUM_LATTICE_NODES,
        encoding=PositionEncoding.UNARY,
    )


@pytest.fixture(scope="module")
def hp_like_interaction_h() -> LigandInteraction:
    return LigandInteraction.hp_like(ligand_hp_type="H")


@pytest.fixture(scope="module")
def hp_like_interaction_p() -> LigandInteraction:
    return LigandInteraction.hp_like(ligand_hp_type="P")


@pytest.fixture(scope="module")
def custom_interaction() -> LigandInteraction:
    return LigandInteraction.custom(
        energy_map={"H": -1.5, "P": -0.2},
        default_energy=0.0,
    )


def _make_builder(
    protein: Protein,
    hp_interaction: HPInteraction,
    distance_map: DistanceMap,
    contact_map: ContactMap,
) -> HamiltonianBuilder:
    return HamiltonianBuilder(
        protein=protein,
        interaction=hp_interaction,
        distance_map=distance_map,
        contact_map=contact_map,
    )


def _expected_turn_qubits(protein: Protein) -> int:
    return (len(protein.main_chain) - 1) * QUBITS_PER_TURN


# ===========================================================================
# Backward compatibility - no ligand
# ===========================================================================


class TestBackwardCompatibility:
    """Ensure that sum_hamiltonians() without ligand arguments is unchanged."""

    def test_no_args_returns_sparse_pauli_op(
        self, protein, hp_interaction, distance_map, contact_map
    ):
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        h = builder.sum_hamiltonians()
        assert isinstance(h, SparsePauliOp)

    def test_no_args_qubit_count_unchanged(
        self, protein, hp_interaction, distance_map, contact_map
    ):
        """Without ligand, qubit count == (N-1)*QUBITS_PER_TURN (or backbone size)."""
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        h = builder.sum_hamiltonians()
        # The qubit count should be at least the protein-only register
        protein_qubits = _expected_turn_qubits(protein)
        assert h.num_qubits is not None
        assert h.num_qubits >= protein_qubits

    def test_explicit_none_identical_to_no_args(
        self, protein, hp_interaction, distance_map, contact_map
    ):
        """sum_hamiltonians(ligand=None, ligand_interaction=None) == sum_hamiltonians()."""
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        h1 = builder.sum_hamiltonians()
        h2 = builder.sum_hamiltonians(ligand=None, ligand_interaction=None)
        diff = (h1 - h2).simplify()
        for coeff in diff.coeffs:
            assert abs(coeff) < 1e-10

    def test_no_ligand_qubit_count_matches_with_none(
        self, protein, hp_interaction, distance_map, contact_map
    ):
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        h_no_arg = builder.sum_hamiltonians()
        h_none = builder.sum_hamiltonians(ligand=None, ligand_interaction=None)
        assert h_no_arg.num_qubits == h_none.num_qubits


# ===========================================================================
# _build_ligand_contact_term - unit tests
# ===========================================================================


class TestBuildLigandContactTerm:
    """Unit tests for _build_ligand_contact_term()."""

    def test_none_returns_zero_operator(
        self, protein, hp_interaction, distance_map, contact_map
    ):
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        h_lig = builder._build_ligand_contact_term(ligand=None, ligand_interaction=None)
        assert isinstance(h_lig, SparsePauliOp)
        simplified = h_lig.simplify()
        for coeff in simplified.coeffs:
            assert abs(coeff) < 1e-10

    def test_none_qubit_count_is_protein_register(
        self, protein, hp_interaction, distance_map, contact_map
    ):
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        h_lig = builder._build_ligand_contact_term(ligand=None, ligand_interaction=None)
        # None ligand → zero op on protein-only register
        # protein_qubits defaults to (N-1)*QUBITS_PER_TURN
        expected = _expected_turn_qubits(protein)
        assert h_lig.num_qubits == expected

    def test_binary_ligand_qubit_count(
        self,
        protein,
        hp_interaction,
        distance_map,
        contact_map,
        ligand_binary,
        hp_like_interaction_h,
    ):
        """With BINARY ligand: n_total = protein_qubits + n_pos."""
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        # Build protein hamiltonians to find protein_qubits
        h_bb = builder._build_backbone_contact_term()
        h_bt = builder._add_backtracking_penalty()
        h_fi = builder._build_external_field_term()
        protein_qubits = max(h_bb.num_qubits, h_bt.num_qubits, h_fi.num_qubits)
        h_lig = builder._build_ligand_contact_term(
            ligand=ligand_binary,
            ligand_interaction=hp_like_interaction_h,
            protein_qubits=protein_qubits,
        )
        expected = protein_qubits + ligand_binary.num_position_qubits
        assert h_lig.num_qubits == expected

    def test_unary_ligand_qubit_count(
        self,
        protein,
        hp_interaction,
        distance_map,
        contact_map,
        ligand_unary,
        hp_like_interaction_h,
    ):
        """With UNARY ligand: n_total = protein_qubits + num_lattice_nodes."""
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        h_bb = builder._build_backbone_contact_term()
        h_bt = builder._add_backtracking_penalty()
        h_fi = builder._build_external_field_term()
        protein_qubits = max(h_bb.num_qubits, h_bt.num_qubits, h_fi.num_qubits)
        h_lig = builder._build_ligand_contact_term(
            ligand=ligand_unary,
            ligand_interaction=hp_like_interaction_h,
            protein_qubits=protein_qubits,
        )
        expected = protein_qubits + ligand_unary.num_position_qubits
        assert h_lig.num_qubits == expected

    def test_returns_sparse_pauli_op(
        self,
        protein,
        hp_interaction,
        distance_map,
        contact_map,
        ligand_binary,
        hp_like_interaction_h,
    ):
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        h_lig = builder._build_ligand_contact_term(
            ligand=ligand_binary, ligand_interaction=hp_like_interaction_h
        )
        assert isinstance(h_lig, SparsePauliOp)

    def test_identity_coefficient_equals_sum_of_energies_hp_hydrophobic(
        self,
        protein,
        hp_interaction,
        distance_map,
        contact_map,
        ligand_binary,
        hp_like_interaction_h,
    ):
        """Variant A: identity coeff of H_ligand = Σ_i E_i (HP, hydrophobic ligand)."""
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        h_lig = builder._build_ligand_contact_term(
            ligand=ligand_binary, ligand_interaction=hp_like_interaction_h
        ).simplify()

        # Compute expected: Σ_i E_ligand(aa_i)
        expected_sum = sum(
            hp_like_interaction_h.get_energy(protein.main_chain.get_symbol_at(i))
            for i in range(len(protein.main_chain))
        )

        # The H_ligand should be a scalar multiple of the identity
        # (Variant A result), so the identity coefficient equals expected_sum.
        identity_coeff = 0.0
        for label, coeff in h_lig.to_list():
            if set(label) == {"I"}:
                identity_coeff += float(coeff.real)
            else:
                assert abs(coeff) < 1e-9, (
                    f"Non-identity term in H_ligand: {label} -> {coeff}"
                )

        assert abs(identity_coeff - expected_sum) < 1e-8, (
            f"Identity coeff {identity_coeff} != expected Σ_i E_i = {expected_sum}"
        )

    def test_identity_coefficient_equals_sum_of_energies_custom(
        self,
        protein,
        hp_interaction,
        distance_map,
        contact_map,
        ligand_binary,
        custom_interaction,
    ):
        """Variant A: identity coeff = Σ_i E_i for custom interaction."""
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        h_lig = builder._build_ligand_contact_term(
            ligand=ligand_binary, ligand_interaction=custom_interaction
        ).simplify()

        expected_sum = sum(
            custom_interaction.get_energy(protein.main_chain.get_symbol_at(i))
            for i in range(len(protein.main_chain))
        )

        identity_coeff = 0.0
        for label, coeff in h_lig.to_list():
            if set(label) == {"I"}:
                identity_coeff += float(coeff.real)
            else:
                assert abs(coeff) < 1e-9, (
                    f"Unexpected non-identity term: {label} -> {coeff}"
                )

        assert abs(identity_coeff - expected_sum) < 1e-8

    def test_polar_ligand_all_zero_energy(
        self,
        protein,
        hp_interaction,
        distance_map,
        contact_map,
        ligand_binary,
        hp_like_interaction_p,
    ):
        """Polar ligand → all energies 0 → H_ligand is effectively zero."""
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        h_lig = builder._build_ligand_contact_term(
            ligand=ligand_binary, ligand_interaction=hp_like_interaction_p
        ).simplify()
        for coeff in h_lig.coeffs:
            assert abs(coeff) < 1e-9

    def test_unary_and_binary_give_same_total_energy(
        self,
        protein,
        hp_interaction,
        distance_map,
        contact_map,
        ligand_binary,
        ligand_unary,
        hp_like_interaction_h,
    ):
        """H_ligand identity coefficient equals Σ_i E_i * norm_factor per encoding.

        For BINARY encoding, Σ_k P_L^{(k)} = I_lig  (coefficient of I = 1).
        For UNARY encoding, Σ_k ½(I - Z_k) has coefficient of I = N/2 in the
        full Hilbert space (not the physical one-hot subspace).  We verify each
        encoding against its own analytically-expected identity coefficient.
        """
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        h_bb = builder._build_backbone_contact_term()
        h_bt = builder._add_backtracking_penalty()
        h_fi = builder._build_external_field_term()
        protein_qubits = max(h_bb.num_qubits, h_bt.num_qubits, h_fi.num_qubits)

        energy_sum = sum(
            hp_like_interaction_h.get_energy(protein.main_chain.get_symbol_at(i))
            for i in range(len(protein.main_chain))
        )

        # BINARY: Σ_k P_L^k = I_lig  → identity coeff = energy_sum * 1
        h_bin = builder._build_ligand_contact_term(
            ligand=ligand_binary,
            ligand_interaction=hp_like_interaction_h,
            protein_qubits=protein_qubits,
        ).simplify()
        bin_sum_coeff = sum(
            float(c.real) for lbl, c in h_bin.to_list() if set(lbl) == {"I"}
        )
        assert abs(bin_sum_coeff - energy_sum) < 1e-8, (
            f"BINARY identity coeff {bin_sum_coeff} != energy_sum {energy_sum}"
        )

        # UNARY: Σ_k ½(I-Z_k) = (N/2)·I in full space  → identity coeff = energy_sum * N/2
        n_nodes = ligand_unary.num_lattice_nodes  # = 4
        unary_norm = n_nodes / 2  # coefficient of I in Σ_k P_L^k
        h_un = builder._build_ligand_contact_term(
            ligand=ligand_unary,
            ligand_interaction=hp_like_interaction_h,
            protein_qubits=protein_qubits,
        ).simplify()
        un_sum_coeff = sum(
            float(c.real) for lbl, c in h_un.to_list() if set(lbl) == {"I"}
        )
        expected_unary_coeff = energy_sum * unary_norm
        assert abs(un_sum_coeff - expected_unary_coeff) < 1e-8, (
            f"UNARY identity coeff {un_sum_coeff} != energy_sum*N/2 = {expected_unary_coeff}"
        )


# ===========================================================================
# sum_hamiltonians with ligand - qubit count tests
# ===========================================================================


class TestSumHamiltoniansWithLigand:
    """Integration tests for sum_hamiltonians() when a ligand is provided."""

    def test_qubit_count_increases_by_position_qubits_binary(
        self,
        protein,
        hp_interaction,
        distance_map,
        contact_map,
        ligand_binary,
        hp_like_interaction_h,
    ):
        """With BINARY ligand, total qubits = protein_qubits + n_pos."""
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        h_no_lig = builder.sum_hamiltonians()
        h_with_lig = builder.sum_hamiltonians(
            ligand=ligand_binary, ligand_interaction=hp_like_interaction_h
        )
        expected_increase = ligand_binary.num_position_qubits
        # protein qubits is max across all protein terms
        protein_qubits = h_no_lig.num_qubits
        assert h_with_lig.num_qubits == protein_qubits + expected_increase

    def test_qubit_count_increases_by_position_qubits_unary(
        self,
        protein,
        hp_interaction,
        distance_map,
        contact_map,
        ligand_unary,
        hp_like_interaction_h,
    ):
        """With UNARY ligand, total qubits = protein_qubits + num_lattice_nodes."""
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        h_no_lig = builder.sum_hamiltonians()
        h_with_lig = builder.sum_hamiltonians(
            ligand=ligand_unary, ligand_interaction=hp_like_interaction_h
        )
        expected_increase = ligand_unary.num_position_qubits  # == num_lattice_nodes
        protein_qubits = h_no_lig.num_qubits
        assert h_with_lig.num_qubits == protein_qubits + expected_increase

    def test_qubit_count_formula_binary(
        self,
        protein,
        hp_interaction,
        distance_map,
        contact_map,
        ligand_binary,
        hp_like_interaction_h,
    ):
        """Explicit formula: n_total >= (N-1)*QUBITS_PER_TURN + n_pos."""
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        h = builder.sum_hamiltonians(
            ligand=ligand_binary, ligand_interaction=hp_like_interaction_h
        )
        min_expected = (
            _expected_turn_qubits(protein) + ligand_binary.num_position_qubits
        )
        assert h.num_qubits is not None
        assert h.num_qubits >= min_expected

    def test_qubit_count_formula_unary(
        self,
        protein,
        hp_interaction,
        distance_map,
        contact_map,
        ligand_unary,
        hp_like_interaction_h,
    ):
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        h = builder.sum_hamiltonians(
            ligand=ligand_unary, ligand_interaction=hp_like_interaction_h
        )
        min_expected = _expected_turn_qubits(protein) + ligand_unary.num_position_qubits
        assert h.num_qubits is not None
        assert h.num_qubits >= min_expected

    def test_returns_sparse_pauli_op(
        self,
        protein,
        hp_interaction,
        distance_map,
        contact_map,
        ligand_binary,
        hp_like_interaction_h,
    ):
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        h = builder.sum_hamiltonians(
            ligand=ligand_binary, ligand_interaction=hp_like_interaction_h
        )
        assert isinstance(h, SparsePauliOp)

    def test_with_custom_interaction(
        self,
        protein,
        hp_interaction,
        distance_map,
        contact_map,
        ligand_binary,
        custom_interaction,
    ):
        """Custom interaction: H_total is a valid SparsePauliOp with correct qubits."""
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        h = builder.sum_hamiltonians(
            ligand=ligand_binary, ligand_interaction=custom_interaction
        )
        assert isinstance(h, SparsePauliOp)
        assert h.num_qubits is not None
        assert h.num_qubits > 0

    def test_ligand_shifts_identity_coefficient(
        self,
        protein,
        hp_interaction,
        distance_map,
        contact_map,
        ligand_binary,
        hp_like_interaction_h,
    ):
        """The ligand term shifts the identity coefficient of H_total by Σ_i E_i."""
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        expected_shift = sum(
            hp_like_interaction_h.get_energy(protein.main_chain.get_symbol_at(i))
            for i in range(len(protein.main_chain))
        )

        h_no_lig = builder.sum_hamiltonians()
        h_with_lig = builder.sum_hamiltonians(
            ligand=ligand_binary, ligand_interaction=hp_like_interaction_h
        )

        # Pad h_no_lig to the larger register for comparison
        h_no_lig_padded = pad_to_n_qubits(h_no_lig, h_with_lig.num_qubits)

        diff = (h_with_lig - h_no_lig_padded).simplify()

        identity_shift = 0.0
        for label, coeff in diff.to_list():
            if set(label) == {"I"}:
                identity_shift += float(coeff.real)
            else:
                assert abs(coeff) < 1e-9, (
                    f"Unexpected non-identity term in diff: {label} -> {coeff}"
                )

        assert abs(identity_shift - expected_shift) < 1e-8, (
            f"Identity shift {identity_shift} != expected Σ_i E_i = {expected_shift}"
        )

    def test_binary_and_unary_same_qubit_increase_different_sizes(
        self,
        protein,
        hp_interaction,
        distance_map,
        contact_map,
        ligand_binary,
        ligand_unary,
        hp_like_interaction_h,
    ):
        """BINARY uses fewer qubits than UNARY for the same number of nodes."""
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        h_no_lig = builder.sum_hamiltonians()
        protein_qubits = h_no_lig.num_qubits
        h_bin = builder.sum_hamiltonians(
            ligand=ligand_binary, ligand_interaction=hp_like_interaction_h
        )
        h_un = builder.sum_hamiltonians(
            ligand=ligand_unary, ligand_interaction=hp_like_interaction_h
        )
        # BINARY: 2 position qubits;  UNARY: 4 position qubits
        assert h_bin.num_qubits == protein_qubits + ligand_binary.num_position_qubits
        assert h_un.num_qubits == protein_qubits + ligand_unary.num_position_qubits
        assert h_bin.num_qubits < h_un.num_qubits


# ===========================================================================
# Error handling
# ===========================================================================


class TestErrorHandling:
    """Verify error cases for sum_hamiltonians and _build_ligand_contact_term."""

    def test_only_ligand_raises_value_error(
        self,
        protein,
        hp_interaction,
        distance_map,
        contact_map,
        ligand_binary,
    ):
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        with pytest.raises(ValueError, match="ligand"):
            builder.sum_hamiltonians(ligand=ligand_binary, ligand_interaction=None)

    def test_only_interaction_raises_value_error(
        self,
        protein,
        hp_interaction,
        distance_map,
        contact_map,
        hp_like_interaction_h,
    ):
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        with pytest.raises(ValueError, match="ligand"):
            builder.sum_hamiltonians(
                ligand=None, ligand_interaction=hp_like_interaction_h
            )
