from collections.abc import Iterator

from logger import get_logger
from protein.bead import Bead
from protein.bead.main_bead import MainBead
from protein.chain import Chain

logger = get_logger()


class MainChain(Chain):
    def __init__(self, protein_sequence: str) -> None:
        super().__init__(protein_sequence=protein_sequence)
        logger.debug(
            "Initializing MainChain with protein sequence: %s", protein_sequence
        )

        self.beads = [
            MainBead(
                symbol=bead,
                index=index,
                parent_chain_len=len(protein_sequence),
            )
            for index, bead in enumerate(protein_sequence)
        ]

    def __iter__(self) -> Iterator[Bead]:
        return iter(self.beads)

    def __getitem__(self, index: int) -> Bead:
        return self.beads[index]

    def __len__(self) -> int:
        return len(self.beads)

    def __str__(self) -> str:
        return "".join(bead.symbol for bead in self.beads)
