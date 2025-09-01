
from protein.bead.main_bead import MainBead
from protein.chain import Chain



class MainChain(Chain):
    def __init__(self, protein_sequence: str) -> None:
        super().__init__(protein_sequence=protein_sequence)

        self.beads = [
            MainBead(symbol=bead, index=i, full_identity=self.full_identity)
            for i, bead in enumerate(protein_sequence)
        ]

    @staticmethod
    def build_turn_qubits() -> None:
        pass
