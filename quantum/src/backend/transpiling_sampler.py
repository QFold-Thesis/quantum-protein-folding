"""Wrapper sampler that handles transpilation for IBM backends."""

from __future__ import annotations

from typing import TYPE_CHECKING

from qiskit import transpile
from qiskit.primitives import BaseSamplerV2

from logger import get_logger

if TYPE_CHECKING:
    from collections.abc import Iterable

    from qiskit.primitives.containers import PrimitiveResult, PubResult, SamplerPubLike
    from qiskit.providers import Backend

logger = get_logger()


class TranspilingSampler(BaseSamplerV2):
    """
    Sampler wrapper that transpiles circuits before execution.

    This is required for IBM Quantum backends since March 2024, as they
    no longer accept non-ISA circuits. The wrapper intercepts circuit
    submissions, transpiles them to the backend's instruction set, and
    passes ISA-compliant circuits to the underlying SamplerV2.

    """

    def __init__(self, sampler: BaseSamplerV2, backend: Backend) -> None:
        """
        Initialize the transpiling sampler.

        Args:
            sampler (BaseSamplerV2): The underlying sampler (e.g., IBM SamplerV2).
            backend (Backend): The backend to transpile circuits for.

        """
        self._sampler = sampler
        self._backend = backend
        self._options = sampler.options
        logger.debug(f"TranspilingSampler initialized for backend: {backend.name}")

    @property
    def options(self):
        """Return the sampler options."""
        return self._options

    def run(
        self, pubs: Iterable[SamplerPubLike], *, shots: int | None = None
    ) -> PrimitiveResult[PubResult]:
        """
        Run the sampler with automatic transpilation.

        Args:
            pubs (Iterable[SamplerPubLike]): An iterable of pub-like objects containing circuits to execute.
            shots (int | None): Number of shots per circuit (optional).

        Returns:
            PrimitiveResult[PubResult]: Results from the underlying sampler with transpiled circuits.

        """
        pub_list: list[SamplerPubLike] = list(pubs)
        transpiled_pubs: list[SamplerPubLike] = []

        for pub in pub_list:
            if hasattr(pub, "circuit"):
                circuit = pub.circuit
            elif isinstance(pub, tuple) and len(pub) > 0:
                circuit = pub[0]
            else:
                circuit = pub

            logger.debug(f"Transpiling circuit with {circuit.num_qubits} qubits")
            transpiled_circuit = transpile(
                circuit,
                backend=self._backend,
                optimization_level=3,
            )
            logger.debug(
                f"Transpiled to {transpiled_circuit.num_qubits} qubits "
                f"({transpiled_circuit.size()} gates)"
            )

            if hasattr(pub, "circuit") or isinstance(pub, tuple):
                transpiled_pubs.append((transpiled_circuit, *pub[1:]))
            else:
                transpiled_pubs.append(transpiled_circuit)

        return self._sampler.run(transpiled_pubs, shots=shots)
