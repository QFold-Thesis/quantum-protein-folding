from constants import EMPTY_SIDECHAIN_PLACEHOLDER
from protein.bead.side_bead import SideBead
from protein.chain import Chain
import logging
from logger import get_logger


logger: logging.Logger = get_logger()


class SideChain(Chain):
    def __init__(self, protein_sequence: str) -> None:
        super().__init__(protein_sequence=protein_sequence)
        logger.debug("Initializing SideChain with protein sequence: %s", protein_sequence)

        self.beads = [
            SideBead(
                symbol=bead,
                index=index,
                parent_chain_len=len(protein_sequence),
            )
            for index, bead in enumerate(protein_sequence)
            if bead != EMPTY_SIDECHAIN_PLACEHOLDER
        ]
