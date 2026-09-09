import numpy as np
import pytest

from analysis import (
    compress_with_layout,
    decode_structure,
    expand_bitstring,
    solve_exactly,
)
from analysis.structure_decoder import (
    GAUGE_FIXED_BITS,
    apply_gauge,
    bits_from_bitstring,
    lattice_basis,
    walk_to_coordinates,
)
from constants import QUBITS_PER_TURN
from enums import TurnDirection
from interaction import LigandInteraction
from particle import Ligand
from tests.helpers import build_hamiltonian

SEQUENCE = "APRLRFY"


def test_bits_are_indexed_by_qubit():
    assert bits_from_bitstring("100") == [0, 0, 1]


def test_lattice_basis_vectors_are_unit_length():
    assert np.allclose(np.linalg.norm(lattice_basis(), axis=1), 1.0)


def test_lattice_basis_vectors_are_mutually_symmetric():
    basis = lattice_basis()
    products = [basis[i] @ basis[j] for i in range(4) for j in range(4) if i != j]

    assert np.allclose(products, products[0])


def test_gauge_pins_the_fixed_qubits():
    bits = apply_gauge([1] * 12)

    for qubit, value in GAUGE_FIXED_BITS.items():
        assert bits[qubit] == value


def test_gauge_leaves_free_qubits_alone():
    bits = apply_gauge([1] * 12)

    assert bits[4] == 1
    assert bits[6] == 1


def test_walk_starts_at_the_origin():
    coordinates = walk_to_coordinates([TurnDirection.DIR_0])

    assert np.allclose(coordinates[0], np.zeros(3))


def test_walk_steps_are_unit_length():
    coordinates = walk_to_coordinates(
        [TurnDirection.DIR_0, TurnDirection.DIR_2, TurnDirection.DIR_1]
    )

    assert np.allclose(np.linalg.norm(np.diff(coordinates, axis=0), axis=1), 1.0)


def test_walk_alternates_the_sublattice_sign():
    forward = walk_to_coordinates([TurnDirection.DIR_0], sign_offset=0)
    flipped = walk_to_coordinates([TurnDirection.DIR_0], sign_offset=1)

    assert np.allclose(forward[-1], -flipped[-1])


def test_decoded_chain_has_one_coordinate_per_residue():
    structure = decode_structure("0" * 12, chain_length=7)

    assert structure.coordinates.shape == (7, 3)
    assert len(structure.turns) == 6


def test_decoded_chain_is_a_lattice_walk():
    structure = decode_structure("0" * 12, chain_length=7)

    steps = np.linalg.norm(np.diff(structure.coordinates, axis=0), axis=1)
    assert np.allclose(steps, 1.0)


def test_no_ligand_means_no_ligand_geometry():
    structure = decode_structure("0" * 12, chain_length=7)

    assert structure.ligand_position is None
    assert structure.claimed_contacts == []
    assert structure.realised_contacts == []


def test_contact_distance_means_true_cartesian_adjacency():
    """The lattice and Cartesian metrics coincide exactly at contact distance.

    Squared distances are held as axis coefficients, which equal the Cartesian
    norm only when the separation is a single lattice step. That is precisely the
    case the contact term scores, so the coupling is geometrically honest.
    """
    ligand = Ligand(num_steps=2)
    builder, hamiltonian = build_hamiltonian(
        SEQUENCE, ligand=ligand, ligand_interaction=LigandInteraction.hp_like("H")
    )
    compressed, kept = compress_with_layout(hamiltonian)
    solution = solve_exactly(compressed)

    structure = decode_structure(
        bitstring=expand_bitstring(solution.bitstring, kept, hamiltonian.num_qubits),
        chain_length=len(SEQUENCE),
        ligand=ligand,
        ligand_walk_offset=builder.ligand_walk_offset,
        ligand_contact_offset=builder.ligand_contact_offset,
    )

    contacts = [
        index
        for index, squared in structure.squared_distances.items()
        if np.isclose(squared, 1.0)
    ]
    assert contacts

    for index in contacts:
        cartesian = np.linalg.norm(
            structure.coordinates[index] - structure.ligand_position
        )
        assert np.isclose(cartesian, 1.0, atol=1e-6)


def test_non_contact_residues_are_further_away_in_cartesian_space():
    ligand = Ligand(num_steps=2)
    builder, hamiltonian = build_hamiltonian(
        SEQUENCE, ligand=ligand, ligand_interaction=LigandInteraction.hp_like("H")
    )
    compressed, kept = compress_with_layout(hamiltonian)
    solution = solve_exactly(compressed)

    structure = decode_structure(
        bitstring=expand_bitstring(solution.bitstring, kept, hamiltonian.num_qubits),
        chain_length=len(SEQUENCE),
        ligand=ligand,
        ligand_walk_offset=builder.ligand_walk_offset,
        ligand_contact_offset=builder.ligand_contact_offset,
    )

    for index, squared in structure.squared_distances.items():
        cartesian = np.linalg.norm(
            structure.coordinates[index] - structure.ligand_position
        )
        if not np.isclose(squared, 1.0):
            assert cartesian > 1.0


def test_realised_contacts_are_a_subset_of_claimed_ones():
    ligand = Ligand(num_steps=2)
    builder, hamiltonian = build_hamiltonian(
        SEQUENCE, ligand=ligand, ligand_interaction=LigandInteraction.hp_like("H")
    )
    compressed, kept = compress_with_layout(hamiltonian)
    solution = solve_exactly(compressed)

    structure = decode_structure(
        bitstring=expand_bitstring(solution.bitstring, kept, hamiltonian.num_qubits),
        chain_length=len(SEQUENCE),
        ligand=ligand,
        ligand_walk_offset=builder.ligand_walk_offset,
        ligand_contact_offset=builder.ligand_contact_offset,
    )

    assert set(structure.realised_contacts) <= set(structure.claimed_contacts)


@pytest.mark.parametrize("chain_length", [5, 6, 7])
def test_turn_count_is_one_less_than_the_chain(chain_length):
    bits = "0" * ((chain_length - 1) * QUBITS_PER_TURN)

    assert len(decode_structure(bits, chain_length).turns) == chain_length - 1
