from constants import EMPTY_SIDECHAIN_PLACEHOLDER
from logger import get_logger
from protein.bead.side_bead import SideBead
from protein.chain import Chain

logger = get_logger()


class SideChain(Chain):
    def __init__(self, protein_sequence: str) -> None:
        super().__init__(protein_sequence=protein_sequence)
        logger.debug(
            "Initializing SideChain with protein sequence: %s", protein_sequence
        )

        self.beads = [
            SideBead(
                symbol=bead,
                index=index,
                parent_chain_len=len(protein_sequence),
            )
            for index, bead in enumerate(protein_sequence)
            if bead != EMPTY_SIDECHAIN_PLACEHOLDER
        ]

    def __str__(self) -> str:
        return "".join(bead.symbol for bead in self.beads)
