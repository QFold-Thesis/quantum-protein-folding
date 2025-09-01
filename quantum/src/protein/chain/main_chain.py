
from collections.abc import Iterator
from protein.bead.main_bead import MainBead
from protein.chain import Chain
from protein.bead import Bead


class MainChain(Chain):
    def __init__(self, protein_sequence: str) -> None:
        super().__init__(protein_sequence=protein_sequence)

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
