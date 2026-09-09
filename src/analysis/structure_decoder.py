"""Decoding of register states into lattice geometry.

Turns a bitstring over the full qubit register into something physical: the
sequence of backbone turns, the Cartesian coordinates of every residue, the
ligand's position, and which residue the ligand claims to bind.

The decoder reads the register with the same conventions the operators are built
with, which lets a decoded structure be scored classically and checked against
the Hamiltonian's own eigenvalue. That cross-check is what makes the ligand
coupling verifiable rather than merely plausible.

Gauge-fixed turns
-----------------
:func:`~utils.qubit_utils.fix_qubits` removes the global rotation freedom by
substituting fixed values for qubits 0, 1, 2, 3 and 5. In operator terms it
replaces ``Z_0, Z_2, Z_3`` with ``+1`` and ``Z_1, Z_5`` with ``-1``; since a turn
qubit is ``(I - Z)/2``, that pins those register bits to 0, 1, 0, 0 and 1
respectively. The decoder applies the same substitution so the geometry it
reports is the one the Hamiltonian actually scored.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from constants import (
    FCC_BASIS,
    MAIN_CHAIN_FIFTH_FIXED_POSITION,
    MAIN_CHAIN_FIXED_POSITIONS,
    QUBITS_PER_TURN,
    SIGN_FLIP_SECOND_QUBIT_INDEX,
    SIGN_FLIP_SIXTH_QUBIT_INDEX,
)
from enums import TurnDirection
from logger import get_logger

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from particle.ligand import Ligand

logger = get_logger()

GAUGE_FIXED_BITS: dict[int, int] = {
    qubit: int(qubit in {SIGN_FLIP_SECOND_QUBIT_INDEX, SIGN_FLIP_SIXTH_QUBIT_INDEX})
    for qubit in [*MAIN_CHAIN_FIXED_POSITIONS, MAIN_CHAIN_FIFTH_FIXED_POSITION]
}


@dataclass(frozen=True)
class DecodedStructure:
    """A conformation read out of a register state.

    Attributes:
        turns (list[TurnDirection]): Backbone turn directions, one per bond.
        coordinates (NDArray[np.float64]): Cartesian residue positions, shape
            ``(chain_length, 3)``.
        ligand_position (NDArray[np.float64] | None): Cartesian ligand position,
            or None when no ligand is present.
        ligand_steps (list[TurnDirection]): Lattice steps taken by the ligand.
        claimed_contacts (list[int]): Residue indices the ligand's contact qubits
            flag as bound.
        squared_distances (dict[int, float]): Squared lattice distance from the
            ligand to each residue it could bind.

    """

    turns: list[TurnDirection]
    coordinates: NDArray[np.float64]
    ligand_position: NDArray[np.float64] | None
    ligand_steps: list[TurnDirection]
    claimed_contacts: list[int]
    squared_distances: dict[int, float]

    @property
    def realised_contacts(self) -> list[int]:
        """list[int]: Residues the ligand both claims and actually touches."""
        return [
            index
            for index in self.claimed_contacts
            if np.isclose(self.squared_distances.get(index, np.inf), 1.0)
        ]


def lattice_basis() -> NDArray[np.float64]:
    """Return the normalised tetrahedral basis vectors.

    Returns:
        NDArray[np.float64]: Four unit vectors spanning the diamond lattice,
        shape ``(4, 3)``.

    """
    basis: NDArray[np.float64] = FCC_BASIS.copy()
    return (basis / np.linalg.norm(basis[0])).astype(np.float64)


def bits_from_bitstring(bitstring: str) -> list[int]:
    """Convert a Qiskit-ordered bitstring into qubit-indexed bits.

    Args:
        bitstring (str): Measured bitstring whose leftmost character is the
            highest-indexed qubit.

    Returns:
        list[int]: Bit value at each qubit index, ascending.

    """
    return [int(char) for char in reversed(bitstring)]


def apply_gauge(bits: list[int]) -> list[int]:
    """Overwrite the gauge-fixed turn qubits with the values the operators assume.

    Args:
        bits (list[int]): Qubit-indexed bits.

    Returns:
        list[int]: A copy with the fixed qubits pinned to their gauge values.

    """
    fixed: list[int] = list(bits)
    for qubit, value in GAUGE_FIXED_BITS.items():
        if qubit < len(fixed):
            fixed[qubit] = value

    return fixed


def decode_step(bits: list[int], qubit_base: int) -> TurnDirection:
    """Read the tetrahedral direction of one lattice step.

    In the dense encoding the two qubits of a step index the direction directly,
    matching the one-hot indicators the operators are built from.

    Args:
        bits (list[int]): Qubit-indexed bits.
        qubit_base (int): Index of the step's first qubit.

    Returns:
        TurnDirection: The direction that step took.

    """
    if len(TurnDirection) == QUBITS_PER_TURN:
        active: list[int] = [
            axis for axis in range(QUBITS_PER_TURN) if bits[qubit_base + axis]
        ]
        return TurnDirection(active[0] if active else 0)

    low: int = bits[qubit_base]
    high: int = bits[qubit_base + 1]
    return TurnDirection(2 * low + high)


def walk_to_coordinates(
    directions: list[TurnDirection], sign_offset: int = 0
) -> NDArray[np.float64]:
    """Trace a lattice walk into Cartesian space.

    Args:
        directions (list[TurnDirection]): Direction of each step.
        sign_offset (int, optional): Shifts the alternating sublattice sign so a
            walk can start on either sublattice. Defaults to 0.

    Returns:
        NDArray[np.float64]: Positions visited, including the origin, of shape
        ``(len(directions) + 1, 3)``.

    """
    basis: NDArray[np.float64] = lattice_basis()
    positions: list[NDArray[np.float64]] = [np.zeros(3)]

    for step, direction in enumerate(directions):
        sublattice_sign: int = (-1) ** (step + sign_offset)
        positions.append(positions[-1] + sublattice_sign * basis[direction.value])

    return np.array(positions)


def decode_structure(
    bitstring: str,
    chain_length: int,
    ligand: Ligand | None = None,
    ligand_walk_offset: int = 0,
    ligand_contact_offset: int = 0,
) -> DecodedStructure:
    """Decode a full-register bitstring into lattice geometry.

    Args:
        bitstring (str): State of the full register in Qiskit ordering.
        chain_length (int): Number of residues in the main chain.
        ligand (Ligand | None, optional): Ligand whose registers should also be
            decoded. Defaults to None.
        ligand_walk_offset (int, optional): First qubit of the ligand's walk
            slice. Defaults to 0.
        ligand_contact_offset (int, optional): First qubit of the ligand's
            contact slice. Defaults to 0.

    Returns:
        DecodedStructure: Turns, coordinates and ligand placement.

    """
    bits: list[int] = apply_gauge(bits_from_bitstring(bitstring))

    turns: list[TurnDirection] = [
        decode_step(bits, QUBITS_PER_TURN * turn) for turn in range(chain_length - 1)
    ]
    coordinates: NDArray[np.float64] = walk_to_coordinates(turns)

    if ligand is None:
        return DecodedStructure(
            turns=turns,
            coordinates=coordinates,
            ligand_position=None,
            ligand_steps=[],
            claimed_contacts=[],
            squared_distances={},
        )

    ligand_steps: list[TurnDirection] = [
        decode_step(bits, ligand_walk_offset + QUBITS_PER_TURN * step)
        for step in range(ligand.num_steps)
    ]
    ligand_position: NDArray[np.float64] = walk_to_coordinates(ligand_steps)[-1]

    eligible: list[int] = ligand.eligible_bead_indices(chain_length)
    claimed_contacts: list[int] = [
        bead_index
        for slot, bead_index in enumerate(eligible)
        if bits[ligand_contact_offset + slot]
    ]

    squared_distances: dict[int, float] = {
        bead_index: float(_axis_squared_distance(turns, ligand_steps, bead_index))
        for bead_index in eligible
    }

    logger.debug(
        "Decoded structure: turns=%s, ligand steps=%s, claimed contacts=%s",
        [turn.value for turn in turns],
        [step.value for step in ligand_steps],
        claimed_contacts,
    )

    return DecodedStructure(
        turns=turns,
        coordinates=coordinates,
        ligand_position=ligand_position,
        ligand_steps=ligand_steps,
        claimed_contacts=claimed_contacts,
        squared_distances=squared_distances,
    )


def _axis_squared_distance(
    turns: list[TurnDirection],
    ligand_steps: list[TurnDirection],
    bead_index: int,
) -> float:
    """Compute the lattice squared distance between a residue and the ligand.

    Uses the axis-coefficient convention shared with
    :func:`~utils.lattice_utils.build_squared_distance`, so the result matches
    what the Hamiltonian scores.

    Args:
        turns (list[TurnDirection]): Backbone turn directions.
        ligand_steps (list[TurnDirection]): Ligand walk directions.
        bead_index (int): Residue to measure to.

    Returns:
        float: Squared lattice distance, equal to 1 for nearest neighbours.

    """
    coefficients: NDArray[np.float64] = np.zeros(len(TurnDirection))

    for step in range(bead_index):
        coefficients[turns[step].value] += (-1) ** step

    for step, direction in enumerate(ligand_steps):
        coefficients[direction.value] -= (-1) ** step

    return float((coefficients**2).sum())
