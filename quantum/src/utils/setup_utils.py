import time
from datetime import UTC, datetime
from typing import TYPE_CHECKING, Any

import numpy as np
from qiskit.circuit.library import real_amplitudes
from qiskit.primitives import StatevectorSampler as Sampler
from qiskit.quantum_info import SparsePauliOp
from qiskit_algorithms import SamplingMinimumEigensolverResult, SamplingVQE
from qiskit_algorithms.optimizers import COBYLA

from builder import HamiltonianBuilder
from constants import INTERACTION_TYPE, OUTPUT_DATA_DIR
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
    """Setup the protein folding system components."""
    protein = Protein(
        main_protein_sequence=main_chain, side_protein_sequence=side_chain
    )

    if INTERACTION_TYPE == InteractionType.MJ:
        interaction: Interaction = MJInteraction()
    elif INTERACTION_TYPE == InteractionType.HP:
        interaction: Interaction = HPInteraction()
    else:
        raise InvalidInteractionTypeError

    contact_map = ContactMap(protein=protein)
    distance_map = DistanceMap(protein=protein)

    return protein, interaction, contact_map, distance_map


def build_and_compress_hamiltonian(
    protein: Protein,
    interaction: Interaction,
    contact_map: ContactMap,
    distance_map: DistanceMap,
) -> tuple[SparsePauliOp, SparsePauliOp]:
    """Build and compress the Hamiltonian."""
    h_builder = HamiltonianBuilder(
        protein=protein,
        interaction=interaction,
        distance_map=distance_map,
        contact_map=contact_map,
    )

    hamiltonian = h_builder.sum_hamiltonians()
    logger.debug(f"Original hamiltonian qubits: {hamiltonian.num_qubits}")

    compressed_h = remove_unused_qubits(hamiltonian)
    logger.debug(f"Compressed hamiltonian qubits: {compressed_h.num_qubits}")

    return hamiltonian, compressed_h


def setup_vqe_optimization(
    num_qubits: int,
) -> tuple[SamplingVQE, list[int], list[float]]:
    """Setup VQE optimization components."""
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
        counts.append(eval_count)
        values.append(mean)

    vqe = SamplingVQE(
        sampler=Sampler(),
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
    """Run the VQE optimization."""
    logger.debug("Starting VQE optimization")

    start_time: float = time.perf_counter()

    raw_results: SamplingMinimumEigensolverResult = vqe.compute_minimum_eigenvalue(
        hamiltonian
    )
    duration: float = time.perf_counter() - start_time
    minutes, seconds = divmod(duration, 60)

    logger.debug(f"VQE optimization completed in {int(minutes)}m {seconds:.2f}s")
    return raw_results


def setup_result_analysis(
    raw_results: SamplingMinimumEigensolverResult, protein: Protein, vqe_iterations: list[int], vqe_energies: list[float]
) -> tuple[ResultInterpreter, ResultVisualizer]:
    timestamp: str = datetime.now(tz=UTC).strftime("%Y_%m_%d-%H_%M_%S")
    dirpath: Path = (
        OUTPUT_DATA_DIR / f"{timestamp}-{protein.main_chain!s}-{protein.side_chain!s}"
    )
    dirpath.mkdir(parents=True, exist_ok=True)

    result_interpreter: ResultInterpreter = ResultInterpreter(
        dirpath=dirpath,
        raw_vqe_results=raw_results,
        protein=protein,
        vqe_iterations=vqe_iterations,
        vqe_energies=vqe_energies,
    )
    result_visualizer: ResultVisualizer = ResultVisualizer(dirpath=dirpath)

    return result_interpreter, result_visualizer
