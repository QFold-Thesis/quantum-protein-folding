from builder import HamiltonianBuilder
from constants import EMPTY_SIDECHAIN_PLACEHOLDER
from contact import ContactMap
from distance import DistanceMap
from interaction import MJInteraction
from protein import Protein


def main() -> None:
    main_chain = "APRLRFY"
    side_chain = EMPTY_SIDECHAIN_PLACEHOLDER * len(main_chain)

    protein = Protein(
        main_protein_sequence=main_chain, side_protein_sequence=side_chain
    )

    mj_interaction = MJInteraction()

    contact_map = ContactMap(protein=protein)

    distance_map = DistanceMap(protein=protein)

    builder = HamiltonianBuilder(
        protein=protein,
        mj=mj_interaction,
        distance_map=distance_map,
        contact_map=contact_map,
    )
    hamiltonian = builder.sum_hamiltonians()
    print(  # noqa: T201
        40 * "-",
        "\n hamiltonian\n",
        hamiltonian,
    )


if __name__ == "__main__":
    main()
