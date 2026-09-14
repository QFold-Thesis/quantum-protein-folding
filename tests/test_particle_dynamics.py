import numpy as np
import pytest
from qiskit.quantum_info import SparsePauliOp

from analysis import decode_structure
from analysis.exact_solver import compress_jointly, to_diagonal
from analysis.ligand_analysis import sweep_particle_environment
from analysis.particle_dynamics import (
    analyse_environment,
    build_coordination_operator,
    build_position_projectors,
)
from constants import LATTICE_COORDINATION_NUMBER, QUBITS_PER_TURN
from enums import LigandCoupling
from exceptions import InvalidOperatorError
from interaction import LigandInteraction
from particle import Ligand, RadialPotential
from tests.helpers import build_hamiltonian, diagonal_of

SEQUENCE = "APRLRFY"


@pytest.fixture(scope="module")
def radial_system():
    ligand = Ligand(num_steps=2)
    builder, hamiltonian = build_hamiltonian(
        SEQUENCE,
        ligand=ligand,
        ligand_interaction=LigandInteraction.hp_like("H", energy_scale=3.0),
        coupling=LigandCoupling.RADIAL_POTENTIAL,
        radial_potential=RadialPotential.square_well({1.0: 1.0, 2.0: 0.6, 3.0: 0.3}),
    )
    return ligand, builder, hamiltonian


# --- joint compression -----------------------------------------------------


def test_joint_compression_keeps_qubits_any_operator_uses():
    left = SparsePauliOp.from_list([("IIZ", 1.0)])
    right = SparsePauliOp.from_list([("ZII", 1.0)])

    _compressed, kept = compress_jointly([left, right])

    assert kept == [0, 2]


def test_joint_compression_preserves_alignment():
    left = SparsePauliOp.from_list([("IIZ", 1.0)])
    right = SparsePauliOp.from_list([("ZII", 1.0)])

    compressed, _ = compress_jointly([left, right])

    assert len(to_diagonal(compressed[0])) == len(to_diagonal(compressed[1]))


def test_joint_compression_rejects_mismatched_widths():
    with pytest.raises(InvalidOperatorError, match="register width"):
        compress_jointly(
            [
                SparsePauliOp.from_list([("Z", 1.0)]),
                SparsePauliOp.from_list([("ZZ", 1.0)]),
            ]
        )


def test_joint_compression_rejects_an_empty_list():
    with pytest.raises(InvalidOperatorError, match="at least one"):
        compress_jointly([])


# --- observables -----------------------------------------------------------


def test_coordination_operator_is_integer_valued(radial_system):
    ligand, builder, hamiltonian = radial_system

    coordination = build_coordination_operator(
        chain_length=len(SEQUENCE),
        ligand=ligand,
        num_qubits=hamiltonian.num_qubits,
        walk_qubit_offset=builder.ligand_walk_offset,
    )
    values = diagonal_of(compress_jointly([coordination])[0][0])

    assert np.allclose(values, np.round(values))


def test_coordination_never_exceeds_the_lattice_coordination(radial_system):
    ligand, builder, hamiltonian = radial_system

    coordination = build_coordination_operator(
        chain_length=len(SEQUENCE),
        ligand=ligand,
        num_qubits=hamiltonian.num_qubits,
        walk_qubit_offset=builder.ligand_walk_offset,
    )
    values = diagonal_of(compress_jointly([coordination])[0][0])

    assert values.max() <= LATTICE_COORDINATION_NUMBER
    assert values.min() >= 0


def test_position_projectors_cover_every_configuration(radial_system):
    ligand, builder, hamiltonian = radial_system

    projectors = build_position_projectors(
        ligand=ligand,
        num_qubits=hamiltonian.num_qubits,
        walk_qubit_offset=builder.ligand_walk_offset,
    )

    assert len(projectors) == 4**ligand.num_steps


def test_position_projectors_partition_the_states(radial_system):
    ligand, builder, hamiltonian = radial_system

    projectors = build_position_projectors(
        ligand=ligand,
        num_qubits=hamiltonian.num_qubits,
        walk_qubit_offset=builder.ligand_walk_offset,
    )
    compressed, _ = compress_jointly(list(projectors.values()))
    total = sum(diagonal_of(operator) for operator in compressed)

    assert np.allclose(total, 1.0)


# --- environment -----------------------------------------------------------


def test_environment_reports_a_consistent_buried_fraction(radial_system):
    ligand, builder, hamiltonian = radial_system

    environment = analyse_environment(
        hamiltonian=hamiltonian,
        chain_length=len(SEQUENCE),
        ligand=ligand,
        walk_qubit_offset=builder.ligand_walk_offset,
    )

    assert environment.buried_fraction == pytest.approx(
        environment.coordination_number / LATTICE_COORDINATION_NUMBER
    )


def test_environment_position_weights_are_normalised(radial_system):
    ligand, builder, hamiltonian = radial_system

    environment = analyse_environment(
        hamiltonian=hamiltonian,
        chain_length=len(SEQUENCE),
        ligand=ligand,
        walk_qubit_offset=builder.ligand_walk_offset,
    )

    assert sum(environment.position_weights.values()) == pytest.approx(1.0)


def test_environment_entropy_is_within_its_bounds(radial_system):
    ligand, builder, hamiltonian = radial_system

    environment = analyse_environment(
        hamiltonian=hamiltonian,
        chain_length=len(SEQUENCE),
        ligand=ligand,
        walk_qubit_offset=builder.ligand_walk_offset,
    )

    assert 0.0 <= environment.position_entropy <= np.log2(4**ligand.num_steps)


def test_rejects_a_non_positive_temperature(radial_system):
    ligand, builder, hamiltonian = radial_system

    with pytest.raises(ValueError, match="temperature"):
        analyse_environment(
            hamiltonian=hamiltonian,
            chain_length=len(SEQUENCE),
            ligand=ligand,
            walk_qubit_offset=builder.ligand_walk_offset,
            temperature=0.0,
        )


# --- phenomena -------------------------------------------------------------


def test_attraction_encapsulates_and_repulsion_expels():
    records = sweep_particle_environment(
        SEQUENCE,
        energy_scales=[3.0, -4.0],
        potential=RadialPotential.square_well({1.0: 1.0, 2.0: 0.6, 3.0: 0.3}),
    )
    attractive, repulsive = records

    assert attractive["regime"] == "attractive"
    assert repulsive["regime"] == "repulsive"
    assert attractive["coordination_number"] > repulsive["coordination_number"]
    assert repulsive["expelled"]


def test_coordination_grows_with_attraction():
    records = sweep_particle_environment(
        SEQUENCE,
        energy_scales=[0.0, 1.0, 3.0],
        potential=RadialPotential.square_well({1.0: 1.0, 2.0: 0.6, 3.0: 0.3}),
    )
    coordinations = [record["coordination_number"] for record in records]

    assert coordinations == sorted(coordinations)


def test_a_neutral_particle_is_reported_as_neutral():
    records = sweep_particle_environment(SEQUENCE, energy_scales=[0.0])

    assert records[0]["regime"] == "neutral"
    assert records[0]["interaction_energy"] == pytest.approx(0.0)


def test_strong_interaction_localises_the_particle():
    """Either sign pins the particle down; indifference leaves it delocalised."""
    records = sweep_particle_environment(
        SEQUENCE,
        energy_scales=[0.0, 3.0],
        potential=RadialPotential.square_well({1.0: 1.0, 2.0: 0.6, 3.0: 0.3}),
    )
    neutral, attractive = records

    assert attractive["position_entropy"] < neutral["position_entropy"]


def test_sweep_reports_the_touching_residues():
    records = sweep_particle_environment(
        SEQUENCE,
        energy_scales=[3.0],
        potential=RadialPotential.square_well({1.0: 1.0, 2.0: 0.6, 3.0: 0.3}),
    )

    touching = records[0]["touching_residues"].split()
    assert len(touching) == int(records[0]["coordination_number"])


def test_ligand_separations_require_a_ligand():
    structure = decode_structure(
        "0" * (len(SEQUENCE) - 1) * QUBITS_PER_TURN, len(SEQUENCE)
    )

    with pytest.raises(ValueError, match="no ligand"):
        structure.ligand_separations()
