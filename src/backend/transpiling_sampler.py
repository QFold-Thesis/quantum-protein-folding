"""Wrapper sampler that handles transpilation for IBM backends."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, cast

from qiskit import transpile
from qiskit.circuit.quantumcircuit import QuantumCircuit
from qiskit.primitives import BasePrimitiveJob, BaseSamplerV2
from qiskit.primitives.containers.sampler_pub_result import SamplerPubResult
from qiskit.providers.backend import BackendV2

from logger import get_logger

if TYPE_CHECKING:
    from collections.abc import Iterable

    from qiskit.primitives.containers import PrimitiveResult, SamplerPubLike

logger = get_logger()


class TranspilingSampler(BaseSamplerV2):
    """Sampler wrapper that transpiles circuits before execution.

    This is required for IBM Quantum backends since March 2024, as they
    no longer accept non-ISA circuits. The wrapper intercepts circuit
    submissions, transpiles them to the backend's instruction set, and
    passes ISA-compliant circuits to the underlying SamplerV2.

    """

    def __init__(self, sampler: BaseSamplerV2, backend: BackendV2) -> None:
        """Initialize the transpiling sampler.

        Args:
            sampler (BaseSamplerV2): The underlying sampler (e.g., IBM SamplerV2).
            backend (Backend): The backend to transpile circuits for.

        """
        self._sampler: BaseSamplerV2 = sampler
        self._backend: BackendV2 = backend
        logger.debug("TranspilingSampler initialized for backend: %s", backend.name)

    def run(
        self, pubs: Iterable[SamplerPubLike], *, shots: int | None = None
    ) -> BasePrimitiveJob[PrimitiveResult[SamplerPubResult], Any]:
        """Run the sampler with automatic transpilation.

        Args:
            pubs (Iterable[SamplerPubLike]): An iterable of pub-like objects containing circuits to execute.
            shots (int | None): Number of shots per circuit (optional).

        Returns:
            PrimitiveResult[PubResult]: Results from the underlying sampler with transpiled circuits.

        Raises:
            TypeError: If an unsupported pub type is encountered.

        """
        logger.debug("Running sampler with automatic transpilation of circuits")

        pub_list: list[SamplerPubLike] = list(pubs)
        transpiled_pubs: list[SamplerPubLike] = []

        for pub in pub_list:
            if hasattr(pub, "circuit"):
                circuit: QuantumCircuit = pub.circuit
            elif isinstance(pub, tuple) and len(pub) > 0:
                circuit: QuantumCircuit = pub[0]
            else:
                circuit: QuantumCircuit = cast(QuantumCircuit, pub)

            logger.debug("Transpiling circuit with %s qubits", circuit.num_qubits)
            transpiled_circuit: QuantumCircuit = transpile(
                circuit,
                backend=self._backend,
                optimization_level=3,
            )
            logger.debug(
                "Transpiled to %s qubits (%s gates)",
                transpiled_circuit.num_qubits,
                transpiled_circuit.size(),
            )

            if hasattr(pub, "circuit") or isinstance(pub, tuple):
                transpiled_pubs.append((transpiled_circuit, *pub[1:]))
            else:
                transpiled_pubs.append(transpiled_circuit)

        logger.debug("All circuits transpiled. Submitting to underlying sampler.")

        return self._sampler.run(transpiled_pubs, shots=shots)
