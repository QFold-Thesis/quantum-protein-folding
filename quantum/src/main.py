from constants import EMPTY_SIDECHAIN_PLACEHOLDER
from interaction import MJInteraction
from protein import Protein


def main() -> None:
    main_chain = "APRLRFY"
    side_chain = EMPTY_SIDECHAIN_PLACEHOLDER * len(main_chain)

    protein = Protein(main_chain, side_chain)

    for bead in protein.main_chain:
        print(bead.symbol, bead.index, bead.turn_qubits)  # noqa: T201

    MJInteraction(protein_sequence=main_chain)


if __name__ == "__main__":
    main()
