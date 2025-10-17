from abc import ABC, abstractmethod
from collections.abc import Iterator

from protein.bead import Bead


class Chain(ABC):
    @abstractmethod
    def __init__(self, protein_sequence: str) -> None:
        self.beads: list[Bead] = []

    def get_symbol_at(self, index: int) -> str:
        """
        Returns the bead symbol at the given chain index.

        """
        return self.beads[index].symbol

    def __iter__(self) -> Iterator[Bead]:
        return iter(self.beads)

    def __getitem__(self, index: int) -> Bead:
        return self.beads[index]

    def __len__(self) -> int:
        return len(self.beads)

    def __str__(self) -> str:
        return "".join(bead.symbol for bead in self.beads)
