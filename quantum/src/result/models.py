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


@dataclass(slots=True)
class AxesLimits:
    x: tuple[float, float]
    y: tuple[float, float]
    z: tuple[float, float]


@dataclass(slots=True)
class RotGifConfig:
    figsize: tuple[int, int]
    elev: float
    azim_start: float
    frames: int
    fps: int


@dataclass(slots=True)
class PlotScene:
    coords_arr: np.ndarray
    coords: list[BeadPosition]
    lattice_points_arr: np.ndarray
    bead_colors: np.ndarray
    mid: np.ndarray
    max_range: float
