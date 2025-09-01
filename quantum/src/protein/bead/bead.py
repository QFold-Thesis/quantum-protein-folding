from abc import ABC, abstractmethod

from qiskit.quantum_info import (  # pyright: ignore[reportMissingTypeStubs]
    SparsePauliOp,
)

from constants import CONFORMATION_ENCODING, QUBITS_PER_TURN
from enums import ConformationEncoding, SubLattice
from exceptions import ConformationEncodingError
from utils.qubit_utils import build_full_identity, build_turn_qubit


class Bead(ABC):
    def __init__(self, symbol: str, index: int, parent_chain_len: int) -> None:
        self.symbol: str = symbol
        self.index: int = index
        self.turn_qubits: tuple[SparsePauliOp, ...] = ()
        self.sublattice: SubLattice = SubLattice.B if index % 2 == 1 else SubLattice.A

        self._num_turn_qubits: int = (parent_chain_len - 1) * QUBITS_PER_TURN
        self._has_turn_qubits: bool = index != parent_chain_len - 1

        self._full_identity: SparsePauliOp = build_full_identity(
            num_qubits=self._num_turn_qubits
        )

        if self._has_turn_qubits:
            self._initialize_turn_qubits()

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
            return
        raise ConformationEncodingError

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
