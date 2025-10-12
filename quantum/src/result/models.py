from dataclasses import dataclass

import numpy as np


@dataclass
class SparseVQEOutput:
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


@dataclass
class BeadPosition:
    index: int
    symbol: str
    x: float
    y: float
    z: float

    @property
    def position(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)
