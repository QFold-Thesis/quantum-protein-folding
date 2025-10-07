from constants import EMPTY_SIDECHAIN_PLACEHOLDER
from logger import get_logger
from utils.setup_utils import (
    build_and_compress_hamiltonian,
    setup_folding_system,
    setup_vqe_optimization,
)

logger = get_logger()


def main() -> None:
    main_chain: str = "APRLRFY"
    side_chain: str = EMPTY_SIDECHAIN_PLACEHOLDER * len(main_chain)

    protein, interaction, contact_map, distance_map = setup_folding_system(
        main_chain, side_chain
    )

    _, compressed_h = build_and_compress_hamiltonian(
        protein=protein,
        interaction=interaction,
        contact_map=contact_map,
        distance_map=distance_map,
    )

    vqe, _, _ = setup_vqe_optimization(compressed_h)

    logger.debug("Starting VQE optimization...")
    _ = vqe.compute_minimum_eigenvalue(compressed_h)


if __name__ == "__main__":
    main()
