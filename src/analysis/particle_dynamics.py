"""Behaviour of a free particle embedded in the chain's lattice.

The contact variant answers one question - which residue does the particle
bind - and answers it by construction, since it is told to bind exactly one.
Once the particle becomes a full bead scored by a radial potential, the
interesting quantities are emergent instead, and this module measures them.

Encapsulation
    A lattice site has exactly four neighbours, so "how enclosed is the
    particle" has an exact answer: the number of those four sites occupied by
    residues. A coordination number of 4 is a chain closed completely around the
    particle; 0 is a particle floating free. Nothing in the Hamiltonian rewards
    coordination directly, so whatever value the ground state reaches is the
    outcome of a competition between wrapping around the particle and keeping
    the chain's own geometry cheap.

Expulsion
    Give the particle a positive interaction energy and the same machinery runs
    backwards: the chain pays to keep it at arm's length, and how far it gets
    pushed is read from the distance to the nearest residue.

Positional stabilisation
    Because the whole spectrum is available, the particle's position
    distribution can be computed exactly rather than sampled - a Boltzmann
    average over every state. A sharply peaked distribution means one site is
    genuinely stabilised; a flat one means the particle is delocalised no matter
    what the ground state happens to be.

Every observable here is built as an operator and evaluated over the whole
spectrum at once, so a sweep costs a matrix of numbers rather than a loop over
states.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from analysis.exact_solver import compress_jointly, to_diagonal
from constants import (
    DEFAULT_TEMPERATURE,
    ENCAPSULATION_THRESHOLD,
    LATTICE_CONTACT_DISTANCE,
    LATTICE_COORDINATION_NUMBER,
    QUBITS_PER_TURN,
)
from logger import get_logger
from utils.lattice_utils import build_chain_position, build_squared_distance
from utils.qubit_utils import build_identity_op, fix_qubits
from utils.spectral_utils import build_diagonal_function

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from qiskit.quantum_info import SparsePauliOp

    from particle.ligand import Ligand

logger = get_logger()


@dataclass(frozen=True)
class ParticleEnvironment:
    """What the chain looks like from the particle's point of view.

    Attributes:
        coordination_number (float): Residues occupying the particle's four
            neighbour sites, in the ground state.
        buried_fraction (float): Coordination number as a share of the four
            available sites.
        nearest_distance (float): Squared lattice distance to the closest residue.
        mean_coordination (float): Boltzmann-averaged coordination number.
        position_entropy (float): Shannon entropy of the particle's position
            distribution, in bits. Zero means one site is fully stabilised.
        position_weights (dict[int, float]): Boltzmann weight of each particle
            walk configuration, keyed by its integer encoding.

    """

    coordination_number: float
    buried_fraction: float
    nearest_distance: float
    mean_coordination: float
    position_entropy: float
    position_weights: dict[int, float]

    @property
    def is_encapsulated(self) -> bool:
        """bool: Whether the chain occupies at least half the neighbour sites."""
        return self.buried_fraction >= ENCAPSULATION_THRESHOLD

    @property
    def is_expelled(self) -> bool:
        """bool: Whether no residue touches the particle at all."""
        return self.coordination_number == 0


def build_coordination_operator(
    chain_length: int,
    ligand: Ligand,
    num_qubits: int,
    walk_qubit_offset: int,
) -> SparsePauliOp:
    """Build the operator counting residues that touch the particle.

    Its eigenvalue on any basis state is that state's coordination number, so a
    single diagonal evaluation scores the whole spectrum.

    Args:
        chain_length (int): Number of residues in the main chain.
        ligand (Ligand): The particle whose neighbourhood is counted.
        num_qubits (int): Width of the register the operator must span.
        walk_qubit_offset (int): First qubit of the particle's walk slice.

    Returns:
        SparsePauliOp: Diagonal operator holding the coordination number.

    """
    ligand_position: list[SparsePauliOp] = ligand.position_operators(
        num_qubits=num_qubits, walk_qubit_offset=walk_qubit_offset
    )

    coordination: SparsePauliOp = build_identity_op(num_qubits, 0.0)

    for bead_index in range(chain_length):
        squared_distance: SparsePauliOp = fix_qubits(
            build_squared_distance(
                build_chain_position(bead_index, num_qubits), ligand_position
            )
        )
        coordination = coordination + build_diagonal_function(
            squared_distance, _touching_indicator
        )

    logger.debug(
        "Built coordination operator over %d residues into %d terms.",
        chain_length,
        len(coordination.paulis),
    )
    return coordination.simplify()


def build_nearest_distance_operator(
    chain_length: int,
    ligand: Ligand,
    num_qubits: int,
    walk_qubit_offset: int,
) -> list[SparsePauliOp]:
    """Build the per-residue squared-distance operators to the particle.

    Args:
        chain_length (int): Number of residues in the main chain.
        ligand (Ligand): The particle to measure to.
        num_qubits (int): Width of the register the operators must span.
        walk_qubit_offset (int): First qubit of the particle's walk slice.

    Returns:
        list[SparsePauliOp]: One squared-distance operator per residue.

    """
    ligand_position: list[SparsePauliOp] = ligand.position_operators(
        num_qubits=num_qubits, walk_qubit_offset=walk_qubit_offset
    )

    return [
        fix_qubits(
            build_squared_distance(
                build_chain_position(bead_index, num_qubits), ligand_position
            )
        )
        for bead_index in range(chain_length)
    ]


def build_position_projectors(
    ligand: Ligand, num_qubits: int, walk_qubit_offset: int
) -> dict[int, SparsePauliOp]:
    """Build a projector onto each distinct particle walk configuration.

    The particle's position is encoded in only a handful of qubits, so every one
    of its configurations can be projected onto explicitly.

    Args:
        ligand (Ligand): The particle whose configurations are enumerated.
        num_qubits (int): Width of the register the operators must span.
        walk_qubit_offset (int): First qubit of the particle's walk slice.

    Returns:
        dict[int, SparsePauliOp]: Projector per configuration, keyed by the
        integer value of the walk qubits.

    """
    from constants import CONFORMATION_ENCODING  # noqa: PLC0415
    from utils.lattice_utils import build_step_indicators  # noqa: PLC0415

    step_indicators: list[tuple[SparsePauliOp, ...]] = [
        build_step_indicators(
            qubit_base=walk_qubit_offset + QUBITS_PER_TURN * step,
            num_qubits=num_qubits,
            encoding=CONFORMATION_ENCODING,
        )
        for step in range(ligand.num_steps)
    ]

    projectors: dict[int, SparsePauliOp] = {}

    for configuration in range(len(step_indicators[0]) ** ligand.num_steps):
        projector: SparsePauliOp | None = None
        remaining: int = configuration

        for step in range(ligand.num_steps):
            direction: int = remaining % len(step_indicators[step])
            remaining //= len(step_indicators[step])
            indicator: SparsePauliOp = step_indicators[step][direction]
            projector = indicator if projector is None else projector @ indicator

        if projector is not None:
            projectors[configuration] = fix_qubits(projector.simplify())

    return projectors


def analyse_environment(
    hamiltonian: SparsePauliOp,
    chain_length: int,
    ligand: Ligand,
    walk_qubit_offset: int,
    temperature: float = DEFAULT_TEMPERATURE,
) -> ParticleEnvironment:
    """Measure how the chain arranges itself around the particle.

    Args:
        hamiltonian (SparsePauliOp): Full Hamiltonian of the coupled system.
        chain_length (int): Number of residues in the main chain.
        ligand (Ligand): The particle placed on the lattice.
        walk_qubit_offset (int): First qubit of the particle's walk slice.
        temperature (float, optional): Temperature for the Boltzmann averages.
            Defaults to DEFAULT_TEMPERATURE.

    Returns:
        ParticleEnvironment: Ground-state and thermal descriptors of the
        particle's surroundings.

    """
    num_qubits: int = int(hamiltonian.num_qubits)

    coordination: SparsePauliOp = build_coordination_operator(
        chain_length=chain_length,
        ligand=ligand,
        num_qubits=num_qubits,
        walk_qubit_offset=walk_qubit_offset,
    )
    distances: list[SparsePauliOp] = build_nearest_distance_operator(
        chain_length=chain_length,
        ligand=ligand,
        num_qubits=num_qubits,
        walk_qubit_offset=walk_qubit_offset,
    )
    projectors: dict[int, SparsePauliOp] = build_position_projectors(
        ligand=ligand, num_qubits=num_qubits, walk_qubit_offset=walk_qubit_offset
    )

    ordered_positions: list[int] = sorted(projectors)
    bundle: list[SparsePauliOp] = [
        hamiltonian,
        coordination,
        *distances,
        *(projectors[key] for key in ordered_positions),
    ]
    compressed, _kept = compress_jointly(bundle)

    energies: NDArray[np.float64] = to_diagonal(compressed[0])
    coordination_values: NDArray[np.float64] = to_diagonal(compressed[1])
    distance_values: NDArray[np.float64] = np.array(
        [to_diagonal(operator) for operator in compressed[2 : 2 + chain_length]]
    )

    ground_state: int = int(np.argmin(energies))
    weights: NDArray[np.float64] = _boltzmann_weights(energies, temperature)

    position_weights: dict[int, float] = {}
    for offset, key in enumerate(ordered_positions):
        occupancy: NDArray[np.float64] = to_diagonal(
            compressed[2 + chain_length + offset]
        )
        position_weights[key] = float((weights * occupancy).sum())

    environment = ParticleEnvironment(
        coordination_number=float(coordination_values[ground_state]),
        buried_fraction=float(
            coordination_values[ground_state] / LATTICE_COORDINATION_NUMBER
        ),
        nearest_distance=float(distance_values[:, ground_state].min()),
        mean_coordination=float((weights * coordination_values).sum()),
        position_entropy=_shannon_entropy(np.array(list(position_weights.values()))),
        position_weights=position_weights,
    )

    logger.info(
        "Particle environment: coordination %.0f/%d, nearest d2 %.1f, entropy %.3f bits",
        environment.coordination_number,
        LATTICE_COORDINATION_NUMBER,
        environment.nearest_distance,
        environment.position_entropy,
    )
    return environment


def _touching_indicator(
    squared_distances: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Return 1 where the squared distance equals the contact shell.

    Args:
        squared_distances (NDArray[np.float64]): Squared lattice distances.

    Returns:
        NDArray[np.float64]: Indicator values.

    """
    return np.isclose(squared_distances, LATTICE_CONTACT_DISTANCE).astype(np.float64)


def _boltzmann_weights(
    energies: NDArray[np.float64], temperature: float
) -> NDArray[np.float64]:
    """Return normalised Boltzmann weights for a spectrum.

    Args:
        energies (NDArray[np.float64]): Energy of every basis state.
        temperature (float): Temperature. Must be positive.

    Returns:
        NDArray[np.float64]: Weights summing to one.

    Raises:
        ValueError: If *temperature* is not positive.

    """
    if temperature <= 0.0:
        msg: str = f"temperature must be positive, got {temperature!r}"
        raise ValueError(msg)

    shifted: NDArray[np.float64] = -(energies - energies.min()) / temperature
    weights: NDArray[np.float64] = np.exp(shifted)

    return weights / weights.sum()


def _shannon_entropy(weights: NDArray[np.float64]) -> float:
    """Return the Shannon entropy of a weight vector, in bits.

    Args:
        weights (NDArray[np.float64]): Non-negative weights.

    Returns:
        float: Entropy in bits. Zero when all weight sits on one entry.

    """
    total: float = float(weights.sum())
    if total <= 0.0:
        return 0.0

    normalised: NDArray[np.float64] = weights / total
    positive: NDArray[np.float64] = normalised[normalised > 0.0]

    return float(-(positive * np.log2(positive)).sum())
