from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

from constants import CONFORMATION_ENCODING, QUBITS_PER_TURN
from enums import SubLattice
from exceptions import ConformationEncodingError
from logger import get_logger
from utils.qubit_utils import build_full_identity, build_turn_qubit

if TYPE_CHECKING:
    from qiskit.quantum_info import (  # pyright: ignore[reportMissingTypeStubs]
        SparsePauliOp,
    )

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

        if None in (CONFORMATION_ENCODING, QUBITS_PER_TURN):
            raise ConformationEncodingError

        self.turn_qubits = tuple(
            build_turn_qubit(
                num_qubits=self._num_turn_qubits,
                z_index=QUBITS_PER_TURN * self.index + i,
            )
            for i in range(QUBITS_PER_TURN)
        )
        logger.debug(
            f"Initialized {len(self.turn_qubits)} turn qubits for Bead {self.symbol} | {self.index} ({CONFORMATION_ENCODING.name} encoding)."
        )

    def turn_funcs(
        self,
    ) -> None | tuple[SparsePauliOp, SparsePauliOp, SparsePauliOp, SparsePauliOp]:
        if len(self.turn_qubits) == 0 or not self._has_turn_qubits:
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
