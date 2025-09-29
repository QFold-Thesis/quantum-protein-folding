from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from protein.bead import Bead


class Chain(ABC):
    @abstractmethod
    def __init__(self, protein_sequence: str) -> None:
        self.beads: list[Bead] = []

    # for now - let's keep turn qubits only in beads

    def get_symbol_at(self, index: int) -> str:
        """
        Returns the bead symbol at the given chain index.

        """
        return self.beads[index].symbol
