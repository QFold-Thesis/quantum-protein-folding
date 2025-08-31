from qiskit.quantum_info import SparsePauliOp  # pyright: ignore[reportMissingTypeStubs]

from protein.bead import Bead


class SideBead(Bead):
    def __init__(self, symbol: str, index: int, chain_length: int) -> None:
        super().__init__(symbol, index, chain_length, is_turning=False)

    def _initate_turn_qubits(self) -> None:
        pass

    def turn_0(self) -> SparsePauliOp:
        msg = "Side beads do not have turn qubits."
        raise NotImplementedError(msg)

    def turn_1(self) -> SparsePauliOp:
        msg = "Side beads do not have turn qubits."
        raise NotImplementedError(msg)

    def turn_2(self) -> SparsePauliOp:
        msg = "Side beads do not have turn qubits."
        raise NotImplementedError(msg)

    def turn_3(self) -> SparsePauliOp:
        msg = "Side beads do not have turn qubits."
        raise NotImplementedError(msg)
