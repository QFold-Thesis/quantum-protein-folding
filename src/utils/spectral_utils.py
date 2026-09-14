"""Building arbitrary functions of a diagonal operator.

A shell-resolved potential asks for something a polynomial cannot cheaply give:
an energy that takes prescribed values on the first few lattice shells and is
*exactly zero* everywhere beyond. Interpolating that with a polynomial in the
squared distance needs one node per reachable shell - degree 14 or more for
residues at the far end of the chain - which both explodes the operator and
leaves a tail that quietly drags the particle toward the chain's centroid,
manufacturing the very encapsulation the model is supposed to measure.

This module takes the other route. Every operator here is diagonal, so it is
just a list of values indexed by basis state, and *any* function of it can be
applied elementwise and converted straight back to Pauli form. A diagonal
operator over ``k`` qubits expands uniquely over Pauli-Z strings as

    d[x] = sum_z c_z * (-1)^popcount(z & x)

and the coefficients come from a Walsh-Hadamard transform, ``c = H d / 2^k``.
The result is exact by construction and has no notion of degree, so a
finite-range well costs no more than a smooth one.

The transform runs on only the qubits the operator actually acts on, which for
a residue near the chain start is a handful, so the cost stays far below the
full register.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from qiskit.quantum_info import Pauli, SparsePauliOp

from analysis.exact_solver import MAX_EXACT_QUBITS, to_diagonal
from exceptions import InvalidOperatorError
from logger import get_logger
from utils.qubit_utils import build_identity_op, find_unused_qubits

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import NDArray

logger = get_logger()

PAULI_COEFFICIENT_TOLERANCE: float = 1e-9


def walsh_hadamard_transform(values: NDArray[np.float64]) -> NDArray[np.float64]:
    """Apply the unnormalised Walsh-Hadamard transform to a vector.

    Args:
        values (NDArray[np.float64]): Vector whose length is a power of two.

    Returns:
        NDArray[np.float64]: The transformed vector.

    Raises:
        ValueError: If the length of *values* is not a power of two.

    """
    length: int = len(values)
    if length & (length - 1):
        msg: str = f"length must be a power of two, got {length}"
        raise ValueError(msg)

    transformed: NDArray[np.float64] = values.astype(np.float64).copy()
    stride: int = 1

    while stride < length:
        for start in range(0, length, stride * 2):
            lower = transformed[start : start + stride].copy()
            upper = transformed[start + stride : start + stride * 2].copy()
            transformed[start : start + stride] = lower + upper
            transformed[start + stride : start + stride * 2] = lower - upper
        stride *= 2

    return transformed


def _restrict_to_acting_qubits(
    operator: SparsePauliOp, acting_qubits: list[int]
) -> SparsePauliOp:
    """Rewrite an operator over only the qubits it acts on.

    Args:
        operator (SparsePauliOp): Diagonal operator to restrict.
        acting_qubits (list[int]): Indices of the qubits it acts on.

    Returns:
        SparsePauliOp: The same operator expressed on a narrower register.

    """
    zeros: NDArray[np.bool_] = np.zeros(len(acting_qubits), dtype=bool)
    restricted: list[Pauli] = [
        Pauli((pauli.z[acting_qubits], zeros)) for pauli in operator.paulis
    ]

    return SparsePauliOp(restricted, operator.coeffs).simplify()


def build_diagonal_function(
    operator: SparsePauliOp,
    function: Callable[[NDArray[np.float64]], NDArray[np.float64]],
    tolerance: float = PAULI_COEFFICIENT_TOLERANCE,
) -> SparsePauliOp:
    """Build ``function(operator)`` exactly for a diagonal operator.

    The function is applied to the operator's eigenvalues and the result is
    converted back into a Pauli-Z expansion, so any shape is representable -
    including one that is identically zero outside a finite range.

    Args:
        operator (SparsePauliOp): Diagonal operator to transform. Must contain
            Pauli-Z and identity factors only.
        function (Callable): Vectorised map applied to the eigenvalues.
        tolerance (float, optional): Pauli coefficients smaller than this are
            dropped. Defaults to PAULI_COEFFICIENT_TOLERANCE.

    Returns:
        SparsePauliOp: Operator equal to *function* applied to *operator*, on the
        same register width.

    Raises:
        InvalidOperatorError: If the operator has no qubit count, or acts on more
            qubits than can be enumerated.

    """
    if operator.num_qubits is None:
        msg: str = "operator.num_qubits is None, cannot build a function of it."
        raise InvalidOperatorError(msg)

    num_qubits: int = int(operator.num_qubits)
    unused: set[int] = set(find_unused_qubits(operator))
    acting_qubits: list[int] = [
        qubit for qubit in range(num_qubits) if qubit not in unused
    ]

    if not acting_qubits:
        constant: float = float(function(to_diagonal(operator))[0])
        return build_identity_op(num_qubits, constant)

    if len(acting_qubits) > MAX_EXACT_QUBITS:
        msg: str = (
            f"Operator acts on {len(acting_qubits)} qubits, above the "
            f"{MAX_EXACT_QUBITS} qubit limit for spectral construction."
        )
        raise InvalidOperatorError(msg)

    restricted: SparsePauliOp = _restrict_to_acting_qubits(operator, acting_qubits)
    eigenvalues: NDArray[np.float64] = to_diagonal(restricted)
    mapped: NDArray[np.float64] = np.asarray(function(eigenvalues), dtype=np.float64)

    coefficients: NDArray[np.float64] = walsh_hadamard_transform(mapped) / len(mapped)
    significant: NDArray[np.int64] = np.flatnonzero(np.abs(coefficients) > tolerance)

    if len(significant) == 0:
        return build_identity_op(num_qubits, 0.0)

    paulis: list[Pauli] = []
    for mask in significant:
        z_table: NDArray[np.bool_] = np.zeros(num_qubits, dtype=bool)
        for position, qubit in enumerate(acting_qubits):
            if (int(mask) >> position) & 1:
                z_table[qubit] = True
        paulis.append(Pauli((z_table, np.zeros(num_qubits, dtype=bool))))

    logger.debug(
        "Built a diagonal function over %d acting qubits into %d Pauli terms.",
        len(acting_qubits),
        len(paulis),
    )
    return SparsePauliOp(paulis, coefficients[significant]).simplify()
