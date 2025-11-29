from typing import TYPE_CHECKING

from constants import EMPTY_SIDECHAIN_PLACEHOLDER
from logger import get_logger
from utils.setup_utils import (
    build_and_compress_hamiltonian,
    run_vqe_optimization,
    setup_folding_system,
    setup_result_analysis,
    setup_vqe_optimization,
)

if TYPE_CHECKING:
    from qiskit_algorithms import SamplingMinimumEigensolverResult

logger = get_logger()


def main() -> None:
    """
    Executes the full quantum protein folding workflow for a sample chain.

    This includes system setup, hamiltonian construction and compression, VQE
    optimization, and result analysis and visualization in 2D and 3D.

    Note:
        The main chain sequence is hardcoded here for demonstration purposes.

    """
    main_chain: str = "APRLRFY"
    side_chain: str = EMPTY_SIDECHAIN_PLACEHOLDER * len(main_chain)

    protein, interaction, contact_map, distance_map = setup_folding_system(
        main_chain=main_chain, side_chain=side_chain
    )

    _, compressed_h = build_and_compress_hamiltonian(
        protein=protein,
        interaction=interaction,
        contact_map=contact_map,
        distance_map=distance_map,
    )

    vqe, counts, values = setup_vqe_optimization(num_qubits=compressed_h.num_qubits)

    raw_results: SamplingMinimumEigensolverResult = run_vqe_optimization(
        vqe=vqe, hamiltonian=compressed_h
    )

    result_interpreter, result_visualizer = setup_result_analysis(
        raw_results=raw_results,
        protein=protein,
        vqe_iterations=counts,
        vqe_energies=values,
    )

    result_interpreter.dump_results_to_files()

    result_visualizer.visualize_3d()
    result_visualizer.visualize_2d()
    result_visualizer.generate_3d_gif()


if __name__ == "__main__":
    main()
