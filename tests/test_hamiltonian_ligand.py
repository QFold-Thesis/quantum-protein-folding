import numpy as np
import pytest

from analysis import (
    compress_with_layout,
    decode_structure,
    expand_bitstring,
    solve_exactly,
)
from builder import HamiltonianBuilder
from constants import QUBITS_PER_TURN
from interaction import LigandInteraction
from particle import Ligand
from tests.helpers import build_hamiltonian, build_system

SEQUENCE = "APRLRFY"


@pytest.fixture(scope="module")
def hydrophobic_ligand():
    return LigandInteraction.hp_like("H")


def solve_with_ligand(sequence, ligand, ligand_interaction):
    builder, hamiltonian = build_hamiltonian(
        sequence, ligand=ligand, ligand_interaction=ligand_interaction
    )
    compressed, kept = compress_with_layout(hamiltonian)
    solution = solve_exactly(compressed)
    structure = decode_structure(
        bitstring=expand_bitstring(solution.bitstring, kept, hamiltonian.num_qubits),
        chain_length=len(sequence),
        ligand=ligand,
        ligand_walk_offset=builder.ligand_walk_offset,
        ligand_contact_offset=builder.ligand_contact_offset,
    )
    return builder, solution, structure


def test_omitting_optional_terms_leaves_the_hamiltonian_untouched():
    _, plain = build_hamiltonian(SEQUENCE)

    assert (
        plain.num_qubits
        == pow(len(SEQUENCE) - 1, 2) + (len(SEQUENCE) - 1) * QUBITS_PER_TURN
    )


def test_ligand_requires_its_interaction_model():
    protein, interaction, contact_map, distance_map = build_system(SEQUENCE)

    with pytest.raises(ValueError, match="ligand_interaction is required"):
        HamiltonianBuilder(
            protein=protein,
            interaction=interaction,
            distance_map=distance_map,
            contact_map=contact_map,
            ligand=Ligand(),
        )


def test_register_layout_places_ligand_above_the_protein(hydrophobic_ligand):
    ligand = Ligand(num_steps=2)
    builder, _ = build_hamiltonian(
        SEQUENCE, ligand=ligand, ligand_interaction=hydrophobic_ligand
    )

    assert builder.ligand_walk_offset == builder.num_protein_qubits
    assert (
        builder.ligand_contact_offset
        == builder.num_protein_qubits + ligand.num_walk_qubits
    )
    assert builder.num_qubits == builder.ligand_contact_offset + len(
        ligand.eligible_bead_indices(len(SEQUENCE))
    )


def test_adding_a_ligand_widens_the_register(hydrophobic_ligand):
    _, plain = build_hamiltonian(SEQUENCE)
    _, with_ligand = build_hamiltonian(
        SEQUENCE, ligand=Ligand(num_steps=2), ligand_interaction=hydrophobic_ligand
    )

    assert with_ligand.num_qubits > plain.num_qubits


@pytest.mark.parametrize("num_steps", [2, 3])
def test_ground_state_never_puts_the_ligand_on_a_residue(num_steps, hydrophobic_ligand):
    ligand = Ligand(num_steps=num_steps)
    _, _, structure = solve_with_ligand(SEQUENCE, ligand, hydrophobic_ligand)

    separations = np.linalg.norm(
        structure.coordinates - structure.ligand_position, axis=1
    )

    assert separations.min() > 0.5


@pytest.mark.parametrize("num_steps", [2, 3])
def test_ground_state_claims_exactly_one_contact(num_steps, hydrophobic_ligand):
    ligand = Ligand(num_steps=num_steps)
    _, _, structure = solve_with_ligand(SEQUENCE, ligand, hydrophobic_ligand)

    assert len(structure.claimed_contacts) == 1


@pytest.mark.parametrize("num_steps", [2, 3])
def test_claimed_contact_is_geometrically_realised(num_steps, hydrophobic_ligand):
    ligand = Ligand(num_steps=num_steps)
    _, _, structure = solve_with_ligand(SEQUENCE, ligand, hydrophobic_ligand)

    assert structure.claimed_contacts == structure.realised_contacts


def test_hydrophobic_ligand_binds_a_hydrophobic_residue(hydrophobic_ligand):
    ligand = Ligand(num_steps=2)
    _, _, structure = solve_with_ligand(SEQUENCE, ligand, hydrophobic_ligand)

    bound = structure.realised_contacts
    assert bound
    assert all(hydrophobic_ligand.get_energy(SEQUENCE[i]) < 0 for i in bound)


def test_binding_lowers_the_energy_by_the_contact_energy(hydrophobic_ligand):
    _, plain = build_hamiltonian(SEQUENCE)
    apo = solve_exactly(compress_with_layout(plain)[0])
    _, holo, _ = solve_with_ligand(SEQUENCE, Ligand(num_steps=2), hydrophobic_ligand)

    assert holo.energy < apo.energy


def test_indifferent_ligand_does_not_change_the_energy():
    neutral = LigandInteraction.custom(energy_map={}, default_energy=0.0)
    _, plain = build_hamiltonian(SEQUENCE)
    apo = solve_exactly(compress_with_layout(plain)[0])
    _, holo, _ = solve_with_ligand(SEQUENCE, Ligand(num_steps=2), neutral)

    assert np.isclose(holo.energy, apo.energy)


def test_stronger_affinity_binds_more_tightly():
    weak = LigandInteraction.hp_like("H", energy_scale=1.0)
    strong = LigandInteraction.hp_like("H", energy_scale=3.0)

    _, weak_solution, _ = solve_with_ligand(SEQUENCE, Ligand(num_steps=2), weak)
    _, strong_solution, _ = solve_with_ligand(SEQUENCE, Ligand(num_steps=2), strong)

    assert strong_solution.energy < weak_solution.energy


def test_ligand_term_is_diagonal(hydrophobic_ligand):
    _, hamiltonian = build_hamiltonian(
        SEQUENCE, ligand=Ligand(num_steps=2), ligand_interaction=hydrophobic_ligand
    )

    assert not hamiltonian.paulis.x.any()
