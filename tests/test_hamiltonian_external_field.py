"""Unit tests for the ExternalField ↔ HamiltonianBuilder integration (Etap A2).

Tests verify:
- external_field=None leaves sum_hamiltonians() output identical to the
  field-free baseline (backward compatibility).
- A uniform field shifts the total Hamiltonian energy by the expected constant.
- A non-uniform field applies per-bead contributions correctly.
- _build_external_field_term() returns a valid SparsePauliOp in all cases.
- The qubit count of the total Hamiltonian is unchanged by the field.

All tests use a minimal but real protein (main_chain="HPPH", side_chain="____")
processed through the full pipeline so that ContactMap, DistanceMap, and
HamiltonianBuilder are exercised with genuine operators.
"""

from __future__ import annotations

import pytest
from qiskit.quantum_info import SparsePauliOp

from src.builder.hamiltonian_builder import HamiltonianBuilder
from src.constants import QUBITS_PER_TURN
from src.contact.contact_map import ContactMap
from src.distance.distance_map import DistanceMap
from src.interaction.hp_interaction import HPInteraction
from src.particle.external_field import ExternalField
from src.protein.protein import Protein

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

MAIN_CHAIN = "HPPHH"  # length 5 - minimum to allow any backbone contacts
SIDE_CHAIN = "_____"


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


def _make_builder(
    protein: Protein,
    hp_interaction: HPInteraction,
    distance_map: DistanceMap,
    contact_map: ContactMap,
    external_field: ExternalField | None = None,
) -> HamiltonianBuilder:
    """Helper: construct a HamiltonianBuilder with the given field."""
    return HamiltonianBuilder(
        protein=protein,
        interaction=hp_interaction,
        distance_map=distance_map,
        contact_map=contact_map,
        external_field=external_field,
    )


# ---------------------------------------------------------------------------
# Backward-compatibility: None field must give the same result as before
# ---------------------------------------------------------------------------


class TestNoneFieldBackwardCompatibility:
    def test_none_field_stored(
        self, protein, hp_interaction, distance_map, contact_map
    ):
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        assert builder.external_field is None

    def test_sum_hamiltonians_identical_without_field(
        self, protein, hp_interaction, distance_map, contact_map
    ):
        """sum_hamiltonians() with external_field=None must equal the baseline."""
        h_no_field = _make_builder(
            protein, hp_interaction, distance_map, contact_map, external_field=None
        ).sum_hamiltonians()

        h_baseline = _make_builder(
            protein, hp_interaction, distance_map, contact_map
        ).sum_hamiltonians()

        # Both SparsePauliOps must represent the same operator.
        diff = (h_no_field - h_baseline).simplify()
        # All coefficients of the difference must be effectively zero.
        for coeff in diff.coeffs:
            assert abs(coeff) < 1e-10, f"Unexpected non-zero coefficient: {coeff}"

    def test_qubit_count_unchanged_by_none_field(
        self, protein, hp_interaction, distance_map, contact_map
    ):
        h_baseline = _make_builder(
            protein, hp_interaction, distance_map, contact_map
        ).sum_hamiltonians()

        h_with_none = _make_builder(
            protein, hp_interaction, distance_map, contact_map, external_field=None
        ).sum_hamiltonians()

        assert h_baseline.num_qubits == h_with_none.num_qubits


# ---------------------------------------------------------------------------
# _build_external_field_term - unit-level checks
# ---------------------------------------------------------------------------


class TestBuildExternalFieldTerm:
    def test_none_field_returns_zero_identity(
        self, protein, hp_interaction, distance_map, contact_map
    ):
        builder = _make_builder(protein, hp_interaction, distance_map, contact_map)
        h_field = builder._build_external_field_term()

        assert isinstance(h_field, SparsePauliOp)
        # Zero operator: all coefficients must be zero after simplification.
        simplified = h_field.simplify()
        for coeff in simplified.coeffs:
            assert abs(coeff) < 1e-10

    def test_uniform_field_term_is_identity_scaled(
        self, protein, hp_interaction, distance_map, contact_map
    ):
        """A uniform field of strength s should produce h_field ≈ N*s · I."""
        strength = -1.5
        chain_len = len(protein.main_chain)
        field = ExternalField.uniform(strength=strength)

        builder = _make_builder(
            protein, hp_interaction, distance_map, contact_map, external_field=field
        )
        h_field = builder._build_external_field_term().simplify()

        assert isinstance(h_field, SparsePauliOp)

        # Expected coefficient: N * strength
        expected_coeff = chain_len * strength
        # All non-identity Pauli terms must be zero; the identity term's
        # coefficient must equal expected_coeff.
        for label, coeff in h_field.to_list():
            if set(label) == {"I"}:
                assert abs(coeff - expected_coeff) < 1e-9, (
                    f"Identity coefficient {coeff} != expected {expected_coeff}"
                )
            else:
                assert abs(coeff) < 1e-10, f"Non-identity term found: {label} {coeff}"

    def test_non_uniform_field_term_coefficients(
        self, protein, hp_interaction, distance_map, contact_map
    ):
        """Non-uniform field: sum of per-bead energies should equal identity coeff."""
        energy_map = {(0,): -2.0, (2,): 1.5}
        chain_len = len(protein.main_chain)
        # Beads not in map default to 0.0
        expected_total = -2.0 + 0.0 + 1.5 + 0.0 * (chain_len - 3)

        field = ExternalField.non_uniform(energy_map, default_energy=0.0)
        builder = _make_builder(
            protein, hp_interaction, distance_map, contact_map, external_field=field
        )
        h_field = builder._build_external_field_term().simplify()

        # Collect the coefficient of the all-I term.
        total_identity_coeff = 0.0
        for label, coeff in h_field.to_list():
            if set(label) == {"I"}:
                total_identity_coeff += float(coeff.real)

        assert abs(total_identity_coeff - expected_total) < 1e-9

    def test_field_term_has_correct_num_qubits(
        self, protein, hp_interaction, distance_map, contact_map
    ):
        """H_field must span the turn-qubit register (n-1)*QUBITS_PER_TURN."""
        expected_qubits = (len(protein.main_chain) - 1) * QUBITS_PER_TURN
        field = ExternalField.uniform(strength=1.0)

        builder = _make_builder(
            protein, hp_interaction, distance_map, contact_map, external_field=field
        )
        h_field = builder._build_external_field_term()

        assert h_field.num_qubits == expected_qubits

    def test_zero_strength_field_is_effectively_zero(
        self, protein, hp_interaction, distance_map, contact_map
    ):
        field = ExternalField.uniform(strength=0.0)
        builder = _make_builder(
            protein, hp_interaction, distance_map, contact_map, external_field=field
        )
        h_field = builder._build_external_field_term().simplify()
        for coeff in h_field.coeffs:
            assert abs(coeff) < 1e-10


# ---------------------------------------------------------------------------
# sum_hamiltonians - integration checks with a real field
# ---------------------------------------------------------------------------


class TestSumHamiltonians:
    def test_uniform_field_shifts_identity_coefficient(
        self, protein, hp_interaction, distance_map, contact_map
    ):
        """Adding a uniform field of strength s adds N·s to the identity coefficient.

        H_with_field - H_no_field must equal N·s · I (at the SparsePauliOp level).
        This is verified by inspecting Pauli coefficients, avoiding O(4^n) matrix
        construction.
        """
        strength = -2.0
        chain_len = len(protein.main_chain)
        expected_shift = chain_len * strength

        h_base = _make_builder(
            protein, hp_interaction, distance_map, contact_map
        ).sum_hamiltonians()

        field = ExternalField.uniform(strength=strength)
        h_with_field = _make_builder(
            protein, hp_interaction, distance_map, contact_map, external_field=field
        ).sum_hamiltonians()

        # Diff must be a pure scalar (all-I term) with coeff = N·s
        diff = (h_with_field - h_base).simplify()
        for label, coeff in diff.to_list():
            if set(label) == {"I"}:
                assert abs(float(coeff.real) - expected_shift) < 1e-8, (
                    f"Identity coeff {coeff} != expected shift {expected_shift}"
                )
            else:
                assert abs(coeff) < 1e-10, (
                    f"Non-identity term in diff: {label} → {coeff}"
                )

    def test_field_does_not_change_qubit_count(
        self, protein, hp_interaction, distance_map, contact_map
    ):
        h_no_field = _make_builder(
            protein, hp_interaction, distance_map, contact_map
        ).sum_hamiltonians()

        field = ExternalField.uniform(strength=5.0)
        h_with_field = _make_builder(
            protein, hp_interaction, distance_map, contact_map, external_field=field
        ).sum_hamiltonians()

        assert h_no_field.num_qubits == h_with_field.num_qubits

    def test_non_uniform_field_same_qubits(
        self, protein, hp_interaction, distance_map, contact_map
    ):
        field = ExternalField.non_uniform({(0,): -3.0, (4,): 1.0})
        h = _make_builder(
            protein, hp_interaction, distance_map, contact_map, external_field=field
        ).sum_hamiltonians()

        assert isinstance(h, SparsePauliOp)
        assert h.num_qubits is not None
        assert h.num_qubits > 0

    def test_non_uniform_field_shifts_identity_by_sum_of_per_bead_energies(
        self, protein, hp_interaction, distance_map, contact_map
    ):
        """Non-uniform field: diff's identity coefficient must equal Σ_i E_field((i,))."""
        chain_len = len(protein.main_chain)
        energy_map = {(0,): -3.0, (2,): 2.5}
        # Beads not in map default to 0.0; compute expected shift manually.
        expected_shift = sum(energy_map.get((i,), 0.0) for i in range(chain_len))

        h_base = _make_builder(
            protein, hp_interaction, distance_map, contact_map
        ).sum_hamiltonians()

        field = ExternalField.non_uniform(energy_map, default_energy=0.0)
        h_with_field = _make_builder(
            protein, hp_interaction, distance_map, contact_map, external_field=field
        ).sum_hamiltonians()

        diff = (h_with_field - h_base).simplify()
        identity_coeff = 0.0
        for label, coeff in diff.to_list():
            if set(label) == {"I"}:
                identity_coeff += float(coeff.real)
            else:
                assert abs(coeff) < 1e-10, (
                    f"Unexpected non-identity term: {label} → {coeff}"
                )

        assert abs(identity_coeff - expected_shift) < 1e-8, (
            f"Identity shift {identity_coeff} != expected {expected_shift}"
        )

    def test_external_field_stored_on_builder(
        self, protein, hp_interaction, distance_map, contact_map
    ):
        field = ExternalField.uniform(strength=1.0)
        builder = _make_builder(
            protein, hp_interaction, distance_map, contact_map, external_field=field
        )
        assert builder.external_field is field
