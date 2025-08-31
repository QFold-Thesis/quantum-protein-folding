from collections.abc import Iterator

from protein.bead import Bead
from protein.bead.main_bead import MainBead
from protein.chain import Chain


class MainChain(Chain):
    def __init__(self, protein_sequence: str) -> None:
        chain_length = len(protein_sequence)
        self.beads = [
            MainBead(
                symbol=symbol,
                index=i,
                chain_length=chain_length,
                is_turning=(i < chain_length - 1),  # last bead is not turning
            )
            for i, symbol in enumerate(protein_sequence)
        ]

    @staticmethod
    def build_turn_qubit() -> None:
        pass

    def __iter__(self) -> Iterator[Bead]:
        return iter(self.beads)

    def __getitem__(self, index: int) -> Bead:
        return self.beads[index]

    def __len__(self):
        return len(self.beads)
