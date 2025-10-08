from dataclasses import dataclass

import numpy as np
from qiskit_algorithms import SamplingMinimumEigensolverResult


@dataclass
class VQEOutput:
    bitstring: str
    probability: float
    state: str
    energy_value: np.complex128

    def __repr__(self) -> str:
        return (
            f"VQEOutput(bitstring={self.bitstring},\n"
            f"  probability={self.probability},\n"
            f"  state={self.state},\n"
            f"  energy_value={self.energy_value})"
        )


def interpret_raw_vqe_output(raw_output: SamplingMinimumEigensolverResult) -> VQEOutput:
    best_measurement = raw_output.best_measurement

    if not best_measurement:
        msg = "No best measurement found in VQE output."
        raise ValueError(msg)

    bitstring: str | None = best_measurement.get("bitstring")
    probability: float | None = best_measurement.get("probability")
    state: str | None = best_measurement.get("state")
    energy_value: np.complex128 | None = best_measurement.get("value")

    if None in (bitstring, probability, state, energy_value):
        msg = "Incomplete best measurement data in VQE output."
        raise ValueError(msg)

    return VQEOutput(
        bitstring=bitstring,
        probability=probability,
        state=state,
        energy_value=energy_value,
    )
