"""
Bead abstraction for quantum protein folding.

Defines the base class `Bead`, representing a single amino acid in the peptide
chain and responsible for constructing turn-related quantum operators used in
the folding model.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

from qiskit.quantum_info import SparsePauliOp

from constants import CONFORMATION_ENCODING, QUBITS_PER_TURN
from enums import SubLattice
from exceptions import ConformationEncodingError
from logger import get_logger
from utils.qubit_utils import build_identity_op, build_turn_qubit

if TYPE_CHECKING:
    from qiskit.quantum_info import SparsePauliOp

logger = get_logger()


class Bead(ABC):
    """
    An abstract class defining a bead of a peptide.

    Attributes:
        symbol (str): One-letter amino acid symbol.
        index (int): Position of the bead in the parent chain.
        turn_qubits (tuple[SparsePauliOp, ...]): Quantum turn qubits associated with this bead.
        sublattice (SubLattice): Sublattice type (A or B) based on bead index.

    """

    def __init__(self, symbol: str, index: int, parent_chain_len: int) -> None:
        """
        Initialize a bead instance.

        Sets symbol, index, sublattice type, and initializes turn qubits if applicable.

        Args:
            symbol (str): One-letter amino acid symbol.
            index (int): Position of the bead in the parent chain.
            parent_chain_len (int): Total number of beads in the parent chain.

        """
        self.symbol: str = symbol
        self.index: int = index
        self.turn_qubits: tuple[SparsePauliOp, ...] = ()
        self.sublattice: SubLattice = SubLattice.B if index % 2 == 1 else SubLattice.A

        self._num_turn_qubits: int = (parent_chain_len - 1) * QUBITS_PER_TURN
        self._has_turn_qubits: bool = index != (parent_chain_len - 1)

        self._full_identity: SparsePauliOp = build_identity_op(
            num_qubits=self._num_turn_qubits
        )

        if self._has_turn_qubits:
            self._initialize_turn_qubits()
        else:
            logger.info(
                "Initialized 0 turn qubits for Bead %s (index: %d) [skipped: last bead in chain]",
                self.symbol,
                self.index,
            )

    def _initialize_turn_qubits(self) -> None:
        """
        Initializes the quantum turn qubits associated with this bead based on the selected conformation encoding.

        For each bead (except the last in the chain), the function builds a set of Pauli operators representing
        directional turns in the protein structure. The number of initialized qubits depends on whether the
        encoding is dense or sparse.

        Raises:
            ConformationEncodingError: If conformation encoding or qubit count is not set.

        """
        if not self._has_turn_qubits:
            return

        if None in (CONFORMATION_ENCODING, QUBITS_PER_TURN):
            raise ConformationEncodingError

        self.turn_qubits = tuple(
            build_turn_qubit(
                num_qubits=self._num_turn_qubits,
                z_index=QUBITS_PER_TURN * self.index + i,
            )
            for i in range(QUBITS_PER_TURN)
        )
        logger.info(
            "Initialized %d turn qubits for Bead %s (index: %d)",
            len(self.turn_qubits),
            self.symbol,
            self.index,
        )

    def turn_funcs(
        self,
    ) -> None | tuple[SparsePauliOp, SparsePauliOp, SparsePauliOp, SparsePauliOp]:
        """
        Return Pauli operators representing directional turns.

        Returns:
            tuple[SparsePauliOp, SparsePauliOp, SparsePauliOp, SparsePauliOp] | None:
                Tuple of Pauli operators for turns (0-3), or `None` if no turn qubits exist.

        """
        if not self.turn_qubits or not self._has_turn_qubits:
            return None
        return (self.turn_0(), self.turn_1(), self.turn_2(), self.turn_3())

    @abstractmethod
    def turn_0(self) -> SparsePauliOp:
        """Return Pauli operator representing turn in direction 0."""
        pass

    @abstractmethod
    def turn_1(self) -> SparsePauliOp:
        """Return Pauli operator representing turn in direction 1."""
        pass

    @abstractmethod
    def turn_2(self) -> SparsePauliOp:
        """Return Pauli operator representing turn in direction 2."""
        pass

    @abstractmethod
    def turn_3(self) -> SparsePauliOp:
        """Return Pauli operator representing turn in direction 3."""
        pass