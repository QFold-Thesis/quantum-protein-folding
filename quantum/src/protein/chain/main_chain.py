from collections.abc import Iterator

from logger import get_logger
from protein.bead import Bead
from protein.bead.main_bead import MainBead
from protein.chain import Chain, SideChain

logger = get_logger()


class MainChain(Chain):
    def __init__(self, protein_sequence: str, side_protein_sequences: list[str]) -> None:
        super().__init__(protein_sequence=protein_sequence)
        logger.debug(
            "Initializing MainChain with protein sequence: %s", protein_sequence
        )

        self.side_protein_sequences = side_protein_sequences
        self.beads = [
            MainBead(
                symbol=bead,
                index=index,
                parent_chain_len=len(protein_sequence),
                side_chain=self._create_side_chain(index)
            )
            for index, bead in enumerate(protein_sequence)
        ]

    def _create_side_chain(self, index: int) -> SideChain:
        side_chain_sequence = self.side_protein_sequences[index]
        return SideChain(protein_sequence=side_chain_sequence)

    def __iter__(self) -> Iterator[Bead]:
        return iter(self.beads)

    def __getitem__(self, index: int) -> Bead:
        return self.beads[index]

    def __len__(self) -> int:
        return len(self.beads)

    def __str__(self) -> str:
        return "".join(bead.symbol for bead in self.beads)
