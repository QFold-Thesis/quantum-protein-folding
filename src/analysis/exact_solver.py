"""Exact diagonalisation of folding Hamiltonians.

Every operator this project builds is diagonal in the computational basis: the
turn qubits, the contact flags and the ligand walk all enter through Pauli-Z
strings only. A folding Hamiltonian is therefore a classical energy function in
disguise, and for the register sizes involved here its full spectrum can simply
be read off.

That makes an exact reference cheaply available, which is what turns a VQE run
into a measurement rather than a hope: the variational energy can be compared
against the true ground state, and the sampled bitstring against the true
minimiser.

The diagonal is evaluated without ever forming a matrix. A Pauli-Z string
contributes ``coeff * (-1)^popcount(z & state)`` to basis state ``state``, so the
whole spectrum is a handful of vectorised parity computations.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from exceptions import InvalidOperatorError
from logger import get_logger
from utils.qubit_utils import find_unused_qubits

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from qiskit.quantum_info import SparsePauliOp

logger = get_logger()

MAX_EXACT_QUBITS: int = 24


@dataclass(frozen=True)
class ExactSolution:
    """Ground state of a diagonal Hamiltonian found by exhaustive evaluation.

    Attributes:
        energy (float): Lowest eigenvalue.
        bitstring (str): Minimising basis state, in Qiskit's ordering where the
            leftmost character is the highest-indexed qubit.
        degeneracy (int): Number of basis states sharing the lowest eigenvalue.
        num_qubits (int): Width of the register that was searched.
        spectrum (NDArray[np.float64]): All eigenvalues, indexed by basis state.

    """

    energy: float
    bitstring: str
    degeneracy: int
    num_qubits: int
    spectrum: NDArray[np.float64] = field(repr=False)

    def rank_of(self, bitstring: str) -> int:
        """Return how many distinct energies lie strictly below a given state.

        A rank of 0 means the state is a true ground state.

        Args:
            bitstring (str): Basis state in Qiskit ordering.

        Returns:
            int: Number of distinct lower energies.

        """
        state: int = int(bitstring, 2)
        distinct: NDArray[np.float64] = np.unique(np.round(self.spectrum, decimals=9))
        return int(np.searchsorted(distinct, round(self.spectrum[state], 9)))


def to_diagonal(operator: SparsePauliOp) -> NDArray[np.float64]:
    """Evaluate a diagonal Pauli operator on every computational basis state.

    Args:
        operator (SparsePauliOp): Operator to evaluate. Must contain Pauli-Z and
            identity factors only.

    Returns:
        NDArray[np.float64]: Eigenvalue for each basis state, indexed by the
        integer value of the state.

    Raises:
        InvalidOperatorError: If the operator has no qubit count, is wider than
            MAX_EXACT_QUBITS, or contains an off-diagonal Pauli factor.

    """
    if operator.num_qubits is None:
        msg: str = "operator.num_qubits is None, cannot evaluate its diagonal."
        raise InvalidOperatorError(msg)

    num_qubits: int = int(operator.num_qubits)
    if num_qubits > MAX_EXACT_QUBITS:
        msg: str = (
            f"Refusing to enumerate {2**num_qubits} states for a {num_qubits}-qubit "
            f"operator (limit is {MAX_EXACT_QUBITS} qubits)."
        )
        raise InvalidOperatorError(msg)

    states: NDArray[np.int64] = np.arange(2**num_qubits, dtype=np.int64)
    diagonal: NDArray[np.float64] = np.zeros(2**num_qubits, dtype=np.float64)

    for pauli, coeff in zip(operator.paulis, operator.coeffs, strict=True):
        if pauli.x.any():
            msg: str = (
                "Operator contains an off-diagonal Pauli factor; exact enumeration "
                "only supports Z-type operators."
            )
            raise InvalidOperatorError(msg)

        parity: NDArray[np.int64] = np.zeros(2**num_qubits, dtype=np.int64)
        for qubit, has_z in enumerate(pauli.z):
            if has_z:
                parity ^= (states >> qubit) & 1

        diagonal += float(coeff.real) * (1 - 2 * parity)

    return diagonal


def solve_exactly(operator: SparsePauliOp) -> ExactSolution:
    """Find the exact ground state of a diagonal Hamiltonian.

    Args:
        operator (SparsePauliOp): Diagonal Hamiltonian to minimise.

    Returns:
        ExactSolution: Ground state energy, a minimising bitstring, its
        degeneracy and the full spectrum.

    """
    diagonal: NDArray[np.float64] = to_diagonal(operator)
    num_qubits: int = int(operator.num_qubits)

    minimum: float = float(diagonal.min())
    minimisers: NDArray[np.int64] = np.flatnonzero(
        np.isclose(diagonal, minimum, rtol=0.0, atol=1e-9)
    )
    best_state: int = int(minimisers[0])

    logger.info(
        "Exact ground state over %d qubits: energy=%.6f, degeneracy=%d",
        num_qubits,
        minimum,
        len(minimisers),
    )

    return ExactSolution(
        energy=minimum,
        bitstring=format(best_state, f"0{num_qubits}b"),
        degeneracy=len(minimisers),
        num_qubits=num_qubits,
        spectrum=diagonal,
    )


def compress_with_layout(operator: SparsePauliOp) -> tuple[SparsePauliOp, list[int]]:
    """Drop qubits the operator never acts on, reporting which ones survived.

    Mirrors :func:`~utils.qubit_utils.remove_unused_qubits` but also returns the
    surviving qubit indices, so a bitstring measured on the compressed register
    can be mapped back onto the full one.

    Args:
        operator (SparsePauliOp): Operator to compress.

    Returns:
        tuple[SparsePauliOp, list[int]]: The compressed operator and the original
        indices of the qubits it retains, in ascending order.

    Raises:
        InvalidOperatorError: If the operator has no qubit count.

    """
    if operator.num_qubits is None:
        msg: str = "operator.num_qubits is None, cannot compress operator."
        raise InvalidOperatorError(msg)

    from utils.qubit_utils import remove_unused_qubits  # noqa: PLC0415

    unused: set[int] = set(find_unused_qubits(operator))
    kept: list[int] = [
        qubit for qubit in range(int(operator.num_qubits)) if qubit not in unused
    ]

    logger.debug(
        "Compressing operator from %d to %d qubits (kept %s).",
        operator.num_qubits,
        len(kept),
        kept,
    )
    return remove_unused_qubits(operator), kept


def expand_bitstring(bitstring: str, kept_qubits: list[int], num_qubits: int) -> str:
    """Map a bitstring measured on a compressed register back to the full one.

    Qubits that were dropped during compression do not affect the energy, so they
    are restored as zeros.

    Args:
        bitstring (str): Measured bitstring in Qiskit ordering, leftmost
            character being the highest-indexed qubit of the compressed register.
        kept_qubits (list[int]): Original indices of the retained qubits, as
            returned by :func:`compress_with_layout`.
        num_qubits (int): Width of the full register.

    Returns:
        str: Bitstring over the full register, in Qiskit ordering.

    Raises:
        ValueError: If the bitstring length does not match *kept_qubits*.

    """
    if len(bitstring) != len(kept_qubits):
        msg: str = (
            f"bitstring of length {len(bitstring)} does not match "
            f"{len(kept_qubits)} kept qubits"
        )
        raise ValueError(msg)

    bits: list[str] = ["0"] * num_qubits
    for position, qubit in enumerate(kept_qubits):
        bits[qubit] = bitstring[len(bitstring) - 1 - position]

    return "".join(reversed(bits))
