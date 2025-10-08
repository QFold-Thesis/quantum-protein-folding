from typing import TYPE_CHECKING

from constants import EMPTY_SIDECHAIN_PLACEHOLDER
from logger import get_logger
from utils.setup_utils import (
    build_and_compress_hamiltonian,
    process_results,
    setup_folding_system,
    setup_vqe_optimization,
)

if TYPE_CHECKING:
    from qiskit_algorithms import SamplingMinimumEigensolverResult

logger = get_logger()


def main() -> None:
    main_chain: str = "APRLRFY"
    side_chain: str = EMPTY_SIDECHAIN_PLACEHOLDER * len(main_chain)

    protein, interaction, contact_map, distance_map = setup_folding_system(
        main_chain=main_chain,
        side_chain=side_chain
    )

    _, compressed_h = build_and_compress_hamiltonian(
        protein=protein,
        interaction=interaction,
        contact_map=contact_map,
        distance_map=distance_map,
    )

    vqe, _, _ = setup_vqe_optimization(compressed_h)

    logger.debug("Starting VQE optimization...")

    raw_results: SamplingMinimumEigensolverResult = vqe.compute_minimum_eigenvalue(
        compressed_h
    )

    process_results(raw_results=raw_results, main_chain=main_chain, side_chain=side_chain)


if __name__ == "__main__":
    main()
