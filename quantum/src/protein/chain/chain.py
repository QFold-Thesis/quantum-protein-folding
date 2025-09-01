from abc import ABC, abstractmethod
from qiskit.quantum_info import SparsePauliOp  # pyright: ignore[reportMissingTypeStubs]
from quantum.src.utils.qubit_utils import build_full_identity
from constants import QUBITS_PER_TURN
from utils.qubit_utils import build_full_identity

class Chain(ABC):
    @abstractmethod
    def __init__(self, protein_sequence: str) -> None:
        possible_turns: int = len(protein_sequence) - 1

        self.num_turn_qubits: int = QUBITS_PER_TURN * possible_turns
        self.full_identity: SparsePauliOp = build_full_identity(num_qubits=self.num_turn_qubits)

    @staticmethod
    @abstractmethod
    def build_turn_qubits() -> None:
        pass
