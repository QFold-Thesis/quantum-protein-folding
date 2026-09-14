import numpy as np
import pytest
from qiskit.quantum_info import SparsePauliOp

from analysis import (
    compress_with_layout,
    decode_structure,
    expand_bitstring,
    solve_exactly,
)
from constants import DEFAULT_HARDCORE_ENERGY
from enums import LigandCoupling
from interaction import LigandInteraction
from particle import Ligand, RadialPotential
from tests.helpers import build_hamiltonian, diagonal_of
from utils.spectral_utils import build_diagonal_function, walsh_hadamard_transform

SEQUENCE = "APRLRFY"


def solve_radial(sequence, ligand, ligand_interaction, potential):
    builder, hamiltonian = build_hamiltonian(
        sequence,
        ligand=ligand,
        ligand_interaction=ligand_interaction,
        coupling=LigandCoupling.RADIAL_POTENTIAL,
        radial_potential=potential,
    )
    compressed, kept = compress_with_layout(hamiltonian)
    solution = solve_exactly(compressed)
    structure = decode_structure(
        bitstring=expand_bitstring(solution.bitstring, kept, hamiltonian.num_qubits),
        chain_length=len(sequence),
        ligand=ligand,
        ligand_walk_offset=builder.ligand_walk_offset,
        ligand_contact_offset=builder.contact_flag_offset,
    )
    return builder, solution, structure


# --- spectral construction -------------------------------------------------


def test_walsh_hadamard_is_its_own_inverse_up_to_scale():
    values = np.array([1.0, -2.0, 3.0, 0.5])

    restored = walsh_hadamard_transform(walsh_hadamard_transform(values)) / len(values)

    assert np.allclose(restored, values)


def test_walsh_hadamard_rejects_non_power_of_two():
    with pytest.raises(ValueError, match="power of two"):
        walsh_hadamard_transform(np.ones(3))


def test_diagonal_function_reproduces_the_identity_map():
    operator = SparsePauliOp.from_list([("ZI", 1.0), ("IZ", 2.0)])

    rebuilt = build_diagonal_function(operator, lambda values: values)

    assert np.allclose(diagonal_of(rebuilt), diagonal_of(operator))


def test_diagonal_function_applies_an_arbitrary_shape():
    operator = SparsePauliOp.from_list([("ZI", 1.0), ("IZ", 2.0)])

    squared = build_diagonal_function(operator, lambda values: values**2)

    assert np.allclose(diagonal_of(squared), diagonal_of(operator) ** 2)


def test_diagonal_function_can_be_exactly_zero_outside_a_window():
    operator = SparsePauliOp.from_list([("ZI", 1.0), ("IZ", 2.0)])
    diagonal = diagonal_of(operator)
    target = diagonal.max()

    windowed = build_diagonal_function(
        operator, lambda values: np.isclose(values, target).astype(float)
    )

    assert np.allclose(diagonal_of(windowed), np.isclose(diagonal, target))


def test_diagonal_function_of_a_constant_operator_is_constant():
    operator = SparsePauliOp.from_list([("II", 3.0)])

    doubled = build_diagonal_function(operator, lambda values: 2 * values)

    assert np.allclose(diagonal_of(doubled), 6.0)


# --- potential profile -----------------------------------------------------


def test_contact_only_profile_reaches_the_first_shell():
    potential = RadialPotential.contact_only()

    assert potential.weight(1.0) == pytest.approx(1.0)
    assert potential.weight(2.0) == pytest.approx(0.0)
    assert potential.cutoff == pytest.approx(1.0)


def test_profile_is_exactly_zero_beyond_the_cutoff():
    potential = RadialPotential.square_well({1.0: 1.0, 2.0: 0.5})

    distances = np.array([0.0, 1.0, 2.0, 3.0, 9.0, 32.0])
    energies = potential.evaluate(distances, contact_energy=-1.0)

    assert np.allclose(energies[3:], 0.0)


def test_hardcore_dominates_the_overlap_shell():
    potential = RadialPotential.square_well({1.0: 1.0})

    energies = potential.evaluate(np.array([0.0]), contact_energy=-1.0)

    assert energies[0] == pytest.approx(DEFAULT_HARDCORE_ENERGY)
    assert energies[0] > 0


def test_decaying_profile_weakens_with_distance():
    potential = RadialPotential.decaying(cutoff=4.0)

    weights = [potential.weight(shell) for shell in (1.0, 2.0, 3.0, 4.0)]

    assert weights == sorted(weights, reverse=True)
    assert potential.weight(5.0) == pytest.approx(0.0)


def test_decaying_profile_rejects_a_cutoff_below_contact():
    with pytest.raises(ValueError, match="cutoff"):
        RadialPotential.decaying(cutoff=0.5)


@pytest.mark.parametrize("hardcore", [float("nan"), float("inf")])
def test_non_finite_hardcore_is_rejected(hardcore):
    with pytest.raises(ValueError, match="finite"):
        RadialPotential(hardcore_energy=hardcore)


def test_negative_shell_is_rejected():
    with pytest.raises(ValueError, match="non-negative"):
        RadialPotential(profile={-1.0: 1.0})


# --- coupling behaviour ----------------------------------------------------


def test_radial_coupling_allocates_no_contact_qubits():
    ligand = Ligand(num_steps=2)
    builder, _ = build_hamiltonian(
        SEQUENCE,
        ligand=ligand,
        ligand_interaction=LigandInteraction.hp_like("H"),
        coupling=LigandCoupling.RADIAL_POTENTIAL,
    )

    assert not builder.has_contact_flags
    assert builder.contact_flag_offset is None
    assert builder.num_qubits == builder.num_protein_qubits + ligand.num_walk_qubits


def test_radial_coupling_is_cheaper_than_contact_flags():
    ligand = Ligand(num_steps=2)
    shared = {
        "ligand": ligand,
        "ligand_interaction": LigandInteraction.hp_like("H"),
    }
    flag_builder, _ = build_hamiltonian(
        SEQUENCE, coupling=LigandCoupling.CONTACT_FLAG, **shared
    )
    radial_builder, _ = build_hamiltonian(
        SEQUENCE, coupling=LigandCoupling.RADIAL_POTENTIAL, **shared
    )

    assert radial_builder.num_qubits < flag_builder.num_qubits


def test_contact_only_radial_matches_the_contact_flag_variant():
    """The two variants describe the same physics when their reach agrees."""
    ligand = Ligand(num_steps=2)
    interaction = LigandInteraction.hp_like("H")

    _, flag_hamiltonian = build_hamiltonian(
        SEQUENCE, ligand=ligand, ligand_interaction=interaction
    )
    flag_energy = solve_exactly(compress_with_layout(flag_hamiltonian)[0]).energy

    _, radial_solution, _ = solve_radial(
        SEQUENCE, ligand, interaction, RadialPotential.contact_only()
    )

    assert np.isclose(radial_solution.energy, flag_energy)


def test_wider_well_binds_more_strongly_than_contact_alone():
    ligand = Ligand(num_steps=2)
    interaction = LigandInteraction.hp_like("H")

    _, narrow, _ = solve_radial(
        SEQUENCE, ligand, interaction, RadialPotential.contact_only()
    )
    _, wide, _ = solve_radial(
        SEQUENCE, ligand, interaction, RadialPotential.square_well({1.0: 1.0, 2.0: 0.6})
    )

    assert wide.energy < narrow.energy


def test_particle_never_occupies_a_residue_site():
    ligand = Ligand(num_steps=2)
    _, _, structure = solve_radial(
        SEQUENCE, ligand, LigandInteraction.hp_like("H"), RadialPotential()
    )

    assert structure.ligand_separations().min() > 0.5


def test_attraction_pulls_residues_into_contact():
    ligand = Ligand(num_steps=2)
    _, _, structure = solve_radial(
        SEQUENCE,
        ligand,
        LigandInteraction.hp_like("H", energy_scale=3.0),
        RadialPotential.square_well({1.0: 1.0, 2.0: 0.6, 3.0: 0.3}),
    )

    assert structure.coordination_number >= 1


def test_repulsion_pushes_every_residue_away():
    ligand = Ligand(num_steps=2)
    _, _, structure = solve_radial(
        SEQUENCE,
        ligand,
        LigandInteraction.hp_like("H", energy_scale=-4.0),
        RadialPotential(),
    )

    assert structure.coordination_number == 0
    assert structure.lattice_neighbours == []


def test_stronger_attraction_does_not_reduce_coordination():
    ligand = Ligand(num_steps=2)
    potential = RadialPotential.square_well({1.0: 1.0, 2.0: 0.6, 3.0: 0.3})

    coordinations = []
    for scale in (0.0, 1.0, 3.0):
        _, _, structure = solve_radial(
            SEQUENCE,
            ligand,
            LigandInteraction.hp_like("H", energy_scale=scale),
            potential,
        )
        coordinations.append(structure.coordination_number)

    assert coordinations == sorted(coordinations)


def test_buried_fraction_is_a_share_of_the_neighbour_sites():
    ligand = Ligand(num_steps=2)
    _, _, structure = solve_radial(
        SEQUENCE, ligand, LigandInteraction.hp_like("H"), RadialPotential()
    )

    assert 0.0 <= structure.buried_fraction <= 1.0
    assert structure.buried_fraction == pytest.approx(structure.coordination_number / 4)
