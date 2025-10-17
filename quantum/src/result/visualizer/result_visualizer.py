from pathlib import Path

import numpy as np
from vedo import Sphere, Line, show

from constants import FCC_BASIS
from result.models import BeadPosition
from numpy.typing import NDArray

class ResultVisualizer:
    def __init__(self, dirpath: Path) -> None:
        self.dirpath = dirpath

    def visualize_3d(self, coordinates: list[BeadPosition]) -> None:
        # points = self._generate_tetrahedral_lattice()
        # actors = self._build_lattice_mesh(points)
        # show(actors, axes=1, bg="white")
        
        actors = []
        for bead in coordinates:
            actors.append(Sphere(bead.position, r=0.15, c="skyblue"))
            actors.append(Sphere(bead.position, r=0.01, c="black").legend(bead.symbol))

        # Connect consecutive beads with bonds
        for i in range(len(coordinates)-1):
            p1 = coordinates[i].position
            p2 = coordinates[i+1].position
            actors.append(Line(p1, p2, c="gray"))

        show(actors, axes=1, bg="white")


    @staticmethod
    def _generate_tetrahedral_lattice(n: int = 2, spacing: float = 1.0) -> NDArray[np.float64]:
        """Generate tetrahedral lattice points up to step n."""
        points: list[NDArray[np.float64]] = []
        for i in range(-n, n + 1):
            for j in range(-n, n + 1):
                for k in range(-n, n + 1):
                    pos = i * FCC_BASIS[0] + j * FCC_BASIS[1] + k * FCC_BASIS[2]
                    pos = pos * spacing / np.linalg.norm(FCC_BASIS[0])  # normalize spacing
                    points.append(tuple(pos))
        return np.array(points)
    
    @staticmethod
    def _build_lattice_mesh(points, radius=0.1, cutoff=1.5) -> list:
        actors = []
        # Add spheres
        for p in points:
            actors.append(Sphere(p, r=radius, c="skyblue"))
        
        # Add bonds (lines between close neighbors)
        for i, p1 in enumerate(points):
            for j, p2 in enumerate(points):
                if i < j:
                    dist = np.linalg.norm(p1 - p2)
                    if dist < cutoff:  # near neighbor
                        actors.append(Line(p1, p2, c="gray", alpha=0.6))
        return actors