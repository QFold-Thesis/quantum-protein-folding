"""Lattice walk algebra shared by the peptide backbone and free particles.

Everything placed on the tetrahedral (diamond/FCC) lattice is described as a
*walk*: a position is the signed sum of unit steps, each step selecting one of
four tetrahedral directions. The backbone is the walk defined by the turn
qubits; a free particle such as a ligand is simply another walk living in its
own slice of the qubit register.

For a walk whose step ``s`` is encoded by the qubits starting at ``base_s``:

    r^a = sum_s (-1)^(s + sign_offset) * ind_a(s)

where ``ind_a(s)`` is the one-hot indicator that step ``s`` went in direction
``a``. Positions are kept in *axis-coefficient* space (four components, one per
tetrahedral direction) rather than Cartesian space, and squared distances follow
the same convention as :class:`~distance.distance_map.DistanceMap`:

    d2(u, v) = sum_a (r_u^a - r_v^a)^2

which evaluates to 1 for lattice nearest neighbours. Every operator produced
here is diagonal in the computational basis, so products and powers are exact
and cheap - no Trotterisation or approximation is involved.

Note that this axis-coefficient metric is not the Cartesian norm. Because the
four tetrahedral basis vectors are not orthogonal, the two are related by

    |v|^2 = (4/3) * sum_a c_a^2 - (1/3) * (sum_a c_a)^2

and they agree exactly when the separation is a single lattice step - which is
the case the contact terms score. Away from contact the two metrics diverge, so
the values here should be read as lattice separations rather than distances in
space.

The ``sign_offset`` argument is what couples a free particle to the chain: a
walk of ``m`` steps ends on the same sublattice as backbone bead ``m``, and on
this lattice only sites of *opposite* sublattice can be nearest neighbours.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from constants import (
    DIST_VECTOR_AXES,
    EMPTY_OP_COEFF,
    QUBITS_PER_TURN,
)
from enums import ConformationEncoding
from exceptions import ConformationEncodingError
from utils.qubit_utils import build_identity_op, build_turn_qubit

if TYPE_CHECKING:
    from qiskit.quantum_info import SparsePauliOp


def build_step_indicators(
    qubit_base: int,
    num_qubits: int,
    encoding: ConformationEncoding,
) -> tuple[SparsePauliOp, ...]:
    """Build the four one-hot direction indicators for a single lattice step.

    Mirrors the turn functions of :class:`~protein.bead.Bead`, but takes the
    qubit base index directly so it can describe a step of any walk, not only a
    backbone bead.

    Args:
        qubit_base (int): Absolute index of the first qubit encoding this step.
        num_qubits (int): Total width of the register the operators live on.
        encoding (ConformationEncoding): Turn encoding in use. DENSE packs a
            direction into two qubits, SPARSE uses one qubit per direction.

    Returns:
        tuple[SparsePauliOp, ...]: Four indicator operators, one per tetrahedral
        direction. Exactly one of them evaluates to 1 on any basis state.

    Raises:
        ConformationEncodingError: If the encoding is not DENSE or SPARSE.

    """
    if encoding == ConformationEncoding.SPARSE:
        return tuple(
            build_turn_qubit(z_index=qubit_base + axis, num_qubits=num_qubits)
            for axis in range(DIST_VECTOR_AXES)
        )

    if encoding != ConformationEncoding.DENSE:
        raise ConformationEncodingError

    identity: SparsePauliOp = build_identity_op(num_qubits)
    low: SparsePauliOp = build_turn_qubit(z_index=qubit_base, num_qubits=num_qubits)
    high: SparsePauliOp = build_turn_qubit(
        z_index=qubit_base + 1, num_qubits=num_qubits
    )

    return (
        ((identity - low) @ (identity - high)).simplify(),
        (high @ (high - low)).simplify(),
        (low @ (low - high)).simplify(),
        (low @ high).simplify(),
    )


def build_walk_position(
    qubit_bases: list[int],
    num_qubits: int,
    sign_offset: int = 0,
    encoding: ConformationEncoding | None = None,
) -> list[SparsePauliOp]:
    """Build the axis-coefficient position operators of a lattice walk.

    The walk starts at the lattice origin and takes one step per entry of
    *qubit_bases*, alternating the sublattice sign at every step.

    Args:
        qubit_bases (list[int]): Absolute qubit base index of each step, in walk
            order. An empty list describes a walk that never leaves the origin.
        num_qubits (int): Total width of the register the operators live on.
        sign_offset (int, optional): Shifts the alternating sublattice sign, so a
            walk can be continued from an arbitrary point along another walk.
            Defaults to 0.
        encoding (ConformationEncoding, optional): Turn encoding in use.
            Defaults to the project-wide CONFORMATION_ENCODING.

    Returns:
        list[SparsePauliOp]: Four operators holding the walk's coefficient along
        each tetrahedral direction.

    """
    resolved_encoding: ConformationEncoding = _resolve_encoding(encoding)

    position: list[SparsePauliOp] = [
        build_identity_op(num_qubits, EMPTY_OP_COEFF) for _ in range(DIST_VECTOR_AXES)
    ]

    for step, qubit_base in enumerate(qubit_bases):
        sublattice_sign: int = (-1) ** (step + sign_offset)
        indicators: tuple[SparsePauliOp, ...] = build_step_indicators(
            qubit_base=qubit_base,
            num_qubits=num_qubits,
            encoding=resolved_encoding,
        )

        for axis, indicator in enumerate(indicators):
            position[axis] = position[axis] + sublattice_sign * indicator

    return [component.simplify() for component in position]


def build_chain_position(
    bead_index: int,
    num_qubits: int,
    encoding: ConformationEncoding | None = None,
) -> list[SparsePauliOp]:
    """Build the position operators of a backbone bead relative to the chain start.

    Bead ``i`` sits at the end of the walk formed by the first ``i`` turns, so
    its qubit bases are the standard turn-qubit slots.

    Args:
        bead_index (int): Index of the bead within the main chain. Bead 0 is the
            origin and yields a zero position.
        num_qubits (int): Total width of the register the operators live on.
        encoding (ConformationEncoding, optional): Turn encoding in use.
            Defaults to the project-wide CONFORMATION_ENCODING.

    Returns:
        list[SparsePauliOp]: Four axis-coefficient position operators.

    """
    return build_walk_position(
        qubit_bases=[QUBITS_PER_TURN * turn for turn in range(bead_index)],
        num_qubits=num_qubits,
        encoding=encoding,
    )


def build_squared_distance(
    position_a: list[SparsePauliOp],
    position_b: list[SparsePauliOp],
) -> SparsePauliOp:
    """Build the squared-distance operator between two lattice positions.

    Uses the same axis-coefficient convention as
    :class:`~distance.distance_map.DistanceMap`, so nearest neighbours evaluate
    to 1 and coincident sites to 0.

    Args:
        position_a (list[SparsePauliOp]): Axis-coefficient operators of the first
            position.
        position_b (list[SparsePauliOp]): Axis-coefficient operators of the
            second position.

    Returns:
        SparsePauliOp: Diagonal operator holding the squared lattice distance.

    """
    num_qubits: int = position_a[0].num_qubits
    squared_distance: SparsePauliOp = build_identity_op(num_qubits, EMPTY_OP_COEFF)

    for component_a, component_b in zip(position_a, position_b, strict=True):
        difference: SparsePauliOp = component_a - component_b
        squared_distance = squared_distance + (difference @ difference)

    return squared_distance.simplify()


def _resolve_encoding(encoding: ConformationEncoding | None) -> ConformationEncoding:
    """Fall back to the project-wide encoding when none was supplied.

    Args:
        encoding (ConformationEncoding | None): Explicitly requested encoding.

    Returns:
        ConformationEncoding: The encoding to build operators with.

    """
    if encoding is not None:
        return encoding

    from constants import CONFORMATION_ENCODING  # noqa: PLC0415

    return CONFORMATION_ENCODING
