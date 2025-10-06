from typing import TYPE_CHECKING, Any

import numpy as np
from qiskit.circuit.library import real_amplitudes
from qiskit.primitives import StatevectorSampler as Sampler
from qiskit.quantum_info import SparsePauliOp
from qiskit_algorithms import SamplingVQE
from qiskit_algorithms.optimizers import COBYLA

from builder import HamiltonianBuilder
from constants import EMPTY_SIDECHAIN_PLACEHOLDER, INTERACTION_TYPE
from contact import ContactMap
from distance import DistanceMap
from enums import InteractionType
from exceptions import InvalidInteractionTypeError, InvalidOperatorError
from interaction import HPInteraction, Interaction, MJInteraction
from logger import get_logger
from protein import Protein
from utils.qubit_utils import remove_unused_qubits

if TYPE_CHECKING:
    from qiskit import QuantumCircuit

logger = get_logger()


def __print_best_results(raw_result: dict[str, Any]) -> None:
    """Debug function only - prints the best results from VQE output."""
    logger.info("\n")
    logger.info("Results:")
    logger.info("=" * 60)

    eigenvalue = float(raw_result.eigenvalue)  # type: ignore  # noqa: PGH003
    logger.info(f"Minimum Eigenvalue: {eigenvalue:.6f}")
    logger.info(f"Folding Energy: {eigenvalue:.6f}")

    if raw_result.best_measurement is not None:  # type: ignore  # noqa: PGH003
        logger.info(f"Bitstring: {raw_result.best_measurement['bitstring']}")  # type: ignore  # noqa: PGH003
        logger.info(f"Probability: {raw_result.best_measurement['probability']:.6f}")  # type: ignore  # noqa: PGH003
        logger.info(f"State: {raw_result.best_measurement['state']}")  # type: ignore  # noqa: PGH003
        logger.info(f"Value: {raw_result.best_measurement['value']}")  # type: ignore  # noqa: PGH003

    logger.info("=" * 60)


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
    raw_result = vqe.compute_minimum_eigenvalue(compressed_h)

    __print_best_results(raw_result)  # type: ignore  # noqa: PGH003


if __name__ == "__main__":
    main()
