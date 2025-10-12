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
from exceptions import InvalidInteractionTypeError, InvalidOperatorError
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
    logger.debug("Original hamiltonian qubits: %d", hamiltonian.num_qubits)

    compressed_h = remove_unused_qubits(hamiltonian)
    logger.debug("Compressed hamiltonian qubits: %d", compressed_h.num_qubits)

    return hamiltonian, compressed_h


def setup_vqe_optimization(
    compressed_h: SparsePauliOp,
) -> tuple[SamplingVQE, list[Any], list[Any]]:
    """Setup VQE optimization components."""
    optimizer = COBYLA(maxiter=50)

    if compressed_h.num_qubits is None:
        msg: str = "Hamiltonian number of qubits is None."
        raise InvalidOperatorError(msg)

    ansatz: QuantumCircuit = real_amplitudes(num_qubits=compressed_h.num_qubits, reps=1)

    counts: list[Any] = []
    values: list[Any] = []

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


def setup_result_analysis(
    raw_results: SamplingMinimumEigensolverResult, protein: Protein
) -> tuple[ResultInterpreter, ResultVisualizer]:
    timestamp: str = datetime.now(tz=UTC).strftime("%Y-%m-%d_%H-%M-%S")
    dirpath: Path = (
        OUTPUT_DATA_DIR / f"{timestamp}_{protein.main_chain!s}-{protein.side_chain!s}"
    )
    dirpath.mkdir(parents=True, exist_ok=True)

    result_interpreter: ResultInterpreter = ResultInterpreter(
        dirpath=dirpath, raw_results=raw_results, protein=protein
    )
    result_visualizer: ResultVisualizer = ResultVisualizer(dirpath=dirpath)

    return result_interpreter, result_visualizer
