from __future__ import annotations

from abc import ABC, abstractmethod

from qiskit.quantum_info import (  # pyright: ignore[reportMissingTypeStubs]
    SparsePauliOp,
)

from constants import CONFORMATION_ENCODING, QUBITS_PER_TURN
from enums import ConformationEncoding, SubLattice
from exceptions import ConformationEncodingError
from logger import get_logger
from utils.qubit_utils import build_full_identity, build_turn_qubit

logger = get_logger()


class Bead(ABC):
    def __init__(self, symbol: str, index: int, parent_chain_len: int) -> None:
        self.symbol: str = symbol
        self.index: int = index
        self.turn_qubits: tuple[SparsePauliOp, ...] = ()
        self.sublattice: SubLattice = SubLattice.B if index % 2 == 1 else SubLattice.A

        self._num_turn_qubits: int = (parent_chain_len - 1) * QUBITS_PER_TURN
        self._has_turn_qubits: bool = index != (parent_chain_len - 1)

        self._full_identity: SparsePauliOp = build_full_identity(
            num_qubits=self._num_turn_qubits
        )

        if self._has_turn_qubits:
            self._initialize_turn_qubits()
        else:
            logger.debug(
                f"Bead {self.symbol} | {self.index} is the last in the chain - skipping turn qubit initialization."
            )

    def _initialize_turn_qubits(self) -> None:
        if not self._has_turn_qubits:
            return

        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            self.turn_qubits = (
                build_turn_qubit(
                    num_qubits=self._num_turn_qubits,
                    z_index=QUBITS_PER_TURN * self.index,
                ),
                build_turn_qubit(
                    num_qubits=self._num_turn_qubits,
                    z_index=QUBITS_PER_TURN * self.index + 1,
                ),
            )
            logger.debug(
                f"Initialized {len(self.turn_qubits)} turn qubits for Bead {self.symbol} | {self.index} ({CONFORMATION_ENCODING.name} encoding)."
            )
            return

        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            self.turn_qubits = (
                build_turn_qubit(
                    num_qubits=self._num_turn_qubits,
                    z_index=QUBITS_PER_TURN * self.index,
                ),
                build_turn_qubit(
                    num_qubits=self._num_turn_qubits,
                    z_index=QUBITS_PER_TURN * self.index + 1,
                ),
                build_turn_qubit(
                    num_qubits=self._num_turn_qubits,
                    z_index=QUBITS_PER_TURN * self.index + 2,
                ),
                build_turn_qubit(
                    num_qubits=self._num_turn_qubits,
                    z_index=QUBITS_PER_TURN * self.index + 3,
                ),
            )
            logger.debug(
                f"Initialized {len(self.turn_qubits)} turn qubits for Bead {self.symbol} | {self.index} ({CONFORMATION_ENCODING.name} encoding)."
            )
            return
        raise ConformationEncodingError

    def turn_funcs(self) -> None | tuple[SparsePauliOp, SparsePauliOp, SparsePauliOp, SparsePauliOp]:
        if self.turn_qubits is None:
            return None
        return (self.turn_0(), self.turn_1(), self.turn_2(), self.turn_3())

    @abstractmethod
    def turn_0(self) -> SparsePauliOp:
        pass

    @abstractmethod
    def turn_1(self) -> SparsePauliOp:
        pass

    @abstractmethod
    def turn_2(self) -> SparsePauliOp:
        pass

    @abstractmethod
    def turn_3(self) -> SparsePauliOp:
        pass
