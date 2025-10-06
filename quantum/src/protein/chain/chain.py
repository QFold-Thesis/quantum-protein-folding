from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from protein.bead import Bead
    from collections.abc import Iterator


class Chain(ABC):
    @abstractmethod
    def __init__(self, protein_sequence: str) -> None:
        self.beads: list[Bead] = []

    def get_symbol_at(self, index: int) -> str:
        """
        Returns the bead symbol at the given chain index.

        """
        return self.beads[index].symbol


    def __getitem__(self, index: int) -> Bead:
        return self.beads[index]

    # Allow iteration directly over the chain
    if TYPE_CHECKING:  # pragma: no cover - for type checkers only
        def __iter__(self) -> Iterator[Bead]: ...
    else:
        def __iter__(self):  # runtime simple iterator
            return iter(self.beads)

    def __len__(self) -> int:
        return len(self.beads)

    def __str__(self) -> str:
        return "".join(bead.symbol for bead in self.beads)
