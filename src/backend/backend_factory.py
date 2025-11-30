"""Factory module for creating quantum backend samplers."""

from __future__ import annotations

from typing import TYPE_CHECKING

from backend.transpiling_sampler import TranspilingSampler
from constants import (
    BACKEND_TYPE,
    IBM_QUANTUM_BACKEND_NAME,
    IBM_QUANTUM_SHOTS,
    IBM_QUANTUM_TOKEN,
)
from enums import BackendType
from exceptions import InvalidBackendError
from logger import get_logger

if TYPE_CHECKING:
    from qiskit.primitives import BaseSamplerV2
    from qiskit.providers import Backend

logger = get_logger()


def get_sampler() -> tuple[BaseSamplerV2, Backend | None]:
    """Get the appropriate sampler based on the configured backend type.

    Returns:
        tuple[BaseSamplerV2, Backend | None]: Configured sampler instance and backend (None for local statevector).

    Raises:
        InvalidBackendError: If the backend type is not supported or configuration is invalid.

    """
    if BACKEND_TYPE == BackendType.LOCAL_STATEVECTOR:
        return _get_local_statevector_sampler(), None
    if BACKEND_TYPE == BackendType.IBM_QUANTUM:
        return _get_ibm_quantum_sampler()

    msg: str = f"Unsupported backend type: {BACKEND_TYPE}"
    raise InvalidBackendError(msg)


def _get_local_statevector_sampler() -> BaseSamplerV2:
    """Get a local statevector sampler for ideal simulation.

    Returns:
        BaseSamplerV2: Local statevector sampler instance.

    """
    from qiskit.primitives import StatevectorSampler

    logger.info("Using local StatevectorSampler (ideal simulation)")
    return StatevectorSampler()


def _get_ibm_quantum_sampler() -> tuple[BaseSamplerV2, Backend]:
    """Get a sampler for IBM Quantum hardware with automatic transpilation.

    Requires qiskit-ibm-runtime package and valid IBM Quantum credentials.

    Returns:
        tuple[BaseSamplerV2, Backend]: Transpiling sampler wrapping IBM SamplerV2 and the backend object.

    Raises:
        InvalidBackendError: If IBM runtime package is not installed or credentials are missing.

    """
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2

    token: str | None = IBM_QUANTUM_TOKEN
    backend_name: str | None = IBM_QUANTUM_BACKEND_NAME

    if not token:
        msg: str = (
            "IBM Quantum token not configured. Set IBM_QUANTUM_TOKEN in constants.py "
            "or as environment variable."
        )
        raise InvalidBackendError(msg)

    if not backend_name:
        msg: str = (
            "IBM Quantum backend name not configured. Set IBM_QUANTUM_BACKEND_NAME in constants.py "
            "or as environment variable."
        )
        raise InvalidBackendError(msg)

    logger.info("Connecting to IBM Quantum service...")
    service = QiskitRuntimeService(channel="ibm_quantum_platform", token=token)

    backend = service.backend(backend_name)
    logger.info("Using IBM Quantum backend: %s", backend_name)
    logger.info("Backend status: %s", backend.status())

    ibm_sampler = SamplerV2(mode=backend)
    ibm_sampler.options.default_shots = IBM_QUANTUM_SHOTS

    sampler = TranspilingSampler(sampler=ibm_sampler, backend=backend)

    logger.info("Configured with %s shots", IBM_QUANTUM_SHOTS)
    logger.info("Circuits will be transpiled automatically before execution")

    return sampler, backend
