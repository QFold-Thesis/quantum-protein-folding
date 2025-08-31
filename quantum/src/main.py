from constants import EMPTY_SIDECHAIN_PLACEHOLDER
from interaction import MJInteraction
from protein import Protein


def main() -> None:
    main_chain = "APRLRFY"
    side_chain = EMPTY_SIDECHAIN_PLACEHOLDER * len(main_chain)

    protein = Protein(main_chain, side_chain)

    for bead in protein.main_chain:
        print(40 * "-")  # noqa: T201
        print(bead.symbol, bead.index, bead.turn_qubits)  # noqa: T201

        if not bead.is_turning:
            continue

        for i in range(4):
            turn_op = getattr(bead, f"turn_{i}")()
            print(f"Turn {i} op:", turn_op)  # noqa: T201
        print(40 * "-")  # noqa: T201

    MJInteraction(protein_sequence=main_chain)


if __name__ == "__main__":
    main()
