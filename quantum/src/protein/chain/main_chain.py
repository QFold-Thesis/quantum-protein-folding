from protein.bead.main_bead import MainBead
from protein.chain import Chain


class MainChain(Chain):
    def __init__(self, protein_sequence: str) -> None:
        self.beads = [
            MainBead(symbol=bead, index=i, chain_length=len(protein_sequence)) for i, bead in enumerate(protein_sequence)
        ]

    @staticmethod
    def build_turn_qubit() -> None:
        pass
