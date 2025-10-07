
from typing import TYPE_CHECKING

from constants import EMPTY_SIDECHAIN_PLACEHOLDER
from logger import get_logger
from utils.result_interpretation_utils import VQEOutput, interpret_raw_vqe_output
from utils.setup_utils import (
    build_and_compress_hamiltonian,
    setup_folding_system,
    setup_vqe_optimization,
)
from utils.visualizer_utils import generate_coords_from_bitstring

if TYPE_CHECKING:
    from qiskit_algorithms import SamplingMinimumEigensolverResult

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
    raw_results: SamplingMinimumEigensolverResult = vqe.compute_minimum_eigenvalue(compressed_h)

    interpreted_results: VQEOutput = interpret_raw_vqe_output(raw_results)
    logger.info(f"VQE optimization results: \n{interpreted_results}")

    coords = generate_coords_from_bitstring(bitstring=interpreted_results.bitstring)
    for idx, (residue, coord) in enumerate(zip(main_chain, coords)):
        logger.info(f"Residue {idx} ({residue})\n{coord}")

if __name__ == "__main__":
    main()
