"""
Module providing setup utilities for protein folding quantum simulation, including
Hamiltonian construction, VQE setup, and result processing.
"""

import time
from datetime import datetime
from typing import TYPE_CHECKING, Any

import numpy as np
from qiskit.circuit.library import real_amplitudes
from qiskit.quantum_info import SparsePauliOp
from qiskit_algorithms import SamplingMinimumEigensolverResult, SamplingVQE
from qiskit_algorithms.optimizers import COBYLA

from backend import get_sampler
from builder import HamiltonianBuilder
from constants import DEFAULT_TIMEZONE, INTERACTION_TYPE, RESULTS_DATA_DIRPATH
from contact import ContactMap
from distance import DistanceMap
from enums import InteractionType
from exceptions import InvalidInteractionTypeError
from interaction import HPInteraction, Interaction, MJInteraction
from logger import get_logger
from protein import Protein
from result.interpreter import ResultInterpreter
from result.visualizer import ResultVisualizer
from utils.qubit_utils import remove_unused_qubits

if TYPE_CHECKING:
    from pathlib import Path

    from qiskit import QuantumCircuit

logger = get_logger()


def setup_folding_system(
    main_chain: str, side_chain: str
) -> tuple[Protein, Interaction, ContactMap, DistanceMap]:
    """
    Setup the protein folding system components.

    Args:
        main_chain (str): Main chain protein sequence.
        side_chain (str): Side chain protein sequence.

    Returns:
        tuple[Protein, Interaction, ContactMap, DistanceMap]: The protein, interaction model,
        contact map, and distance map.

    Raises:
        InvalidInteractionTypeError: If the interaction type is invalid (class not inheriting from Interaction).

    """
    if INTERACTION_TYPE == InteractionType.MJ:
        interaction: Interaction = MJInteraction()
    elif INTERACTION_TYPE == InteractionType.HP:
        interaction: Interaction = HPInteraction()
    else:
        raise InvalidInteractionTypeError

    protein = Protein(
        main_protein_sequence=main_chain,
        side_protein_sequence=side_chain,
        valid_symbols=interaction.valid_symbols,
    )

    contact_map = ContactMap(protein=protein)
    distance_map = DistanceMap(protein=protein)

    return protein, interaction, contact_map, distance_map


def build_and_compress_hamiltonian(
    protein: Protein,
    interaction: Interaction,
    contact_map: ContactMap,
    distance_map: DistanceMap,
) -> tuple[SparsePauliOp, SparsePauliOp]:
    """
    Build and compress the final Hamiltonian for the protein folding system.

    Args:
        protein (Protein): The protein instance.
        interaction (Interaction): The interaction model.
        contact_map (ContactMap): The contact map.
        distance_map (DistanceMap): The distance map.

    Returns:
        tuple[SparsePauliOp, SparsePauliOp]: The original and compressed Hamiltonians

    """
    h_builder = HamiltonianBuilder(
        protein=protein,
        interaction=interaction,
        distance_map=distance_map,
        contact_map=contact_map,
    )

    hamiltonian = h_builder.sum_hamiltonians()
    logger.info("Original hamiltonian qubits: %s", hamiltonian.num_qubits)

    compressed_h = remove_unused_qubits(hamiltonian)
    logger.info("Compressed hamiltonian qubits: %s", compressed_h.num_qubits)

    return hamiltonian, compressed_h


def setup_vqe_optimization(
    num_qubits: int,
) -> tuple[SamplingVQE, list[int], list[float]]:
    """
    Setup the VQE optimization process.

    Args:
        num_qubits (int): Number of qubits for the ansatz.

    Returns:
        tuple[SamplingVQE, list[int], list[float]]: The VQE instance, evaluation counts (iterations), and their respective energy values.

    """
    optimizer = COBYLA(maxiter=50)

    ansatz: QuantumCircuit = real_amplitudes(num_qubits=num_qubits, reps=1)

    counts: list[int] = []
    values: list[float] = []

    def _store_intermediate_result(
        eval_count: int,
        _parameters: np.ndarray[Any, Any],
        mean: float,
        _std: dict[str, Any],
    ) -> None:
        """Callback to store intermediate VQE results."""
        counts.append(eval_count)
        values.append(mean)

    sampler, backend = get_sampler()
    if backend is not None:
        logger.debug("Using backend: %s", backend.name)
        logger.debug("Transpilation will be handled by the sampler during execution")

    vqe = SamplingVQE(
        sampler=sampler,
        ansatz=ansatz,
        optimizer=optimizer,
        aggregation=0.1,
        callback=_store_intermediate_result,
    )

    return vqe, counts, values


def run_vqe_optimization(
    vqe: SamplingVQE,
    hamiltonian: SparsePauliOp,
) -> SamplingMinimumEigensolverResult:
    """
    Run the VQE optimization process.

    Args:
        vqe (SamplingVQE): The VQE instance.
        hamiltonian (SparsePauliOp): The Hamiltonian to optimize.

    Returns:
        SamplingMinimumEigensolverResult: The raw results from the VQE optimization.

    """
    logger.debug("Starting VQE optimization")

    start_time: float = time.perf_counter()

    raw_results: SamplingMinimumEigensolverResult = vqe.compute_minimum_eigenvalue(
        hamiltonian
    )
    duration: float = time.perf_counter() - start_time
    minutes, seconds = divmod(duration, 60)

    logger.info("VQE optimization completed in %sm %.2fs", int(minutes), seconds)
    return raw_results


def setup_result_analysis(
    raw_results: SamplingMinimumEigensolverResult,
    protein: Protein,
    vqe_iterations: list[int],
    vqe_energies: list[float],
) -> tuple[ResultInterpreter, ResultVisualizer]:
    """
    Setup the result analysis components.

    Args:
        raw_results (SamplingMinimumEigensolverResult): The raw results from the VQE optimization.
        protein (Protein): The protein instance.
        vqe_iterations (list[int]): The VQE evaluation counts (iterations).
        vqe_energies (list[float]): The VQE energy values.

    Returns:
        tuple[ResultInterpreter, ResultVisualizer]: The result interpreter and visualizer instances.

    """
    timestamp: str = datetime.now(tz=DEFAULT_TIMEZONE).strftime("%Y_%m_%d-%H_%M_%S")
    RESULTS_DATA_DIRPATH.mkdir(parents=True, exist_ok=True)

    dirpath: Path = (
        RESULTS_DATA_DIRPATH / f"{timestamp}-{protein.main_chain!s}-{protein.side_chain!s}"
    )
    dirpath.mkdir(parents=True, exist_ok=True)

    result_interpreter: ResultInterpreter = ResultInterpreter(
        dirpath=dirpath,
        raw_vqe_results=raw_results,
        protein=protein,
        vqe_iterations=vqe_iterations,
        vqe_energies=vqe_energies,
    )

    result_visualizer: ResultVisualizer = ResultVisualizer(
        dirpath=dirpath,
        turn_sequence=result_interpreter.turn_sequence,
        coordinates_3d=result_interpreter.coordinates_3d,
        main_main_contacts_detected=result_interpreter.main_main_contacts_detected,
    )

    return result_interpreter, result_visualizer
