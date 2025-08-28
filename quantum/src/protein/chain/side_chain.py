from constants import EMPTY_SIDECHAIN_PLACEHOLDER
from protein.bead.side_bead import SideBead
from protein.chain import Chain


class SideChain(Chain):
    def __init__(self, protein_sequence: str) -> None:
        super().__init__(protein_sequence=protein_sequence)

        self.beads = [
            SideBead(full_identity=self.full_identity, symbol=bead, index=i)
            for i, bead in enumerate(protein_sequence)
            if bead != EMPTY_SIDECHAIN_PLACEHOLDER
        ]

    @staticmethod
    def build_turn_qubit() -> None:
        _msg: str = "SideChains are not yet implemented!"
        raise NotImplementedError(_msg)
