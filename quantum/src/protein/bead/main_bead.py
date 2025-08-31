from typing import NoReturn

from qiskit.quantum_info import SparsePauliOp  # pyright: ignore[reportMissingTypeStubs]

from constants import CONFORMATION_ENCODING
from enums import ConformationEncoding
from protein.bead import Bead
from utils.qubit_utils import build_full_identity, build_turn_qubit


class MainBead(Bead):
    def __init__(
        self, symbol: str, index: int, chain_length: int, *, is_turning: bool
    ) -> None:
        super().__init__(symbol, index, chain_length, is_turning=is_turning)

        if is_turning:
            self._initate_turn_qubits()
            if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
                # pre-build indicator functions for reuse
                self._turn_fun_0 = self._build_turn_indicator_fun_0()
                self._turn_fun_1 = self._build_turn_indicator_fun_1()
                self._turn_fun_2 = self._build_turn_indicator_fun_2()
                self._turn_fun_3 = self._build_turn_indicator_fun_3()
        else:
            self.turn_qubits: tuple[SparsePauliOp, ...] = ()

    def _initate_turn_qubits(self) -> None:
        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            turn_qubit_0 = build_turn_qubit(4 * (self.chain_length - 1), 4 * self.index)
            turn_qubit_1 = build_turn_qubit(
                4 * (self.chain_length - 1), 4 * self.index + 1
            )
            turn_qubit_2 = build_turn_qubit(
                4 * (self.chain_length - 1), 4 * self.index + 2
            )
            turn_qubit_3 = build_turn_qubit(
                4 * (self.chain_length - 1), 4 * self.index + 3
            )
            self.turn_qubits = (
                turn_qubit_0,
                turn_qubit_1,
                turn_qubit_2,
                turn_qubit_3,
            )

        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            turn_qubit_0 = build_turn_qubit(2 * (self.chain_length - 1), 2 * self.index)
            turn_qubit_1 = build_turn_qubit(
                2 * (self.chain_length - 1), 2 * self.index + 1
            )
            self.turn_qubits = (
                turn_qubit_0,
                turn_qubit_1,
            )

    # ---------------- INDICATOR FUNS (dense only) ----------------
    def _build_turn_indicator_fun_0(self) -> SparsePauliOp:
        full_id = build_full_identity(2 * (self.chain_length - 1))
        return (
            (full_id - self.turn_qubits[0]) @ (full_id - self.turn_qubits[1])
        ).simplify()

    def _build_turn_indicator_fun_1(self) -> SparsePauliOp:
        return (
            self.turn_qubits[1] @ (self.turn_qubits[1] - self.turn_qubits[0])
        ).simplify()

    def _build_turn_indicator_fun_2(self) -> SparsePauliOp:
        return (
            self.turn_qubits[0] @ (self.turn_qubits[0] - self.turn_qubits[1])
        ).simplify()

    def _build_turn_indicator_fun_3(self) -> SparsePauliOp:
        return (self.turn_qubits[0] @ self.turn_qubits[1]).simplify()

    # ---------------- TURN ACCESSORS ----------------

    def _not_implemented(self, turn: int) -> NoReturn:
        msg = f"Invalid encoding for turn_{turn}"
        raise NotImplementedError(msg)

    def turn_0(self) -> SparsePauliOp:
        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            return self.turn_qubits[0]
        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            return self._turn_fun_0
        return self._not_implemented(0)

    def turn_1(self) -> SparsePauliOp:
        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            return self.turn_qubits[1]
        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            return self._turn_fun_1
        return self._not_implemented(1)

    def turn_2(self) -> SparsePauliOp:
        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            return self.turn_qubits[2]
        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            return self._turn_fun_2
        return self._not_implemented(2)

    def turn_3(self) -> SparsePauliOp:
        if CONFORMATION_ENCODING == ConformationEncoding.SPARSE:
            return self.turn_qubits[3]
        if CONFORMATION_ENCODING == ConformationEncoding.DENSE:
            return self._turn_fun_3
        return self._not_implemented(3)
