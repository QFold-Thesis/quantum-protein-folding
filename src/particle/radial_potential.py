"""Distance-resolved potential between a free particle and the chain.

The contact model asks a single yes/no question - is the particle touching this
residue? - and answers it with a flag qubit. A radial potential asks a richer
one: *how far* is the particle from this residue, and what does that cost? The
answer is a profile over lattice shells,

    V(residue, d2) = E(residue) * profile(d2)  +  hardcore * [d2 == 0]

so the interaction matrix gains a second axis. Instead of one number per residue
there is now one number per (residue, shell) pair, which is what lets the model
express a particle that is attracted at touching distance, mildly attracted one
shell out, and forbidden from the residue's own site.

Finite range is the essential property. Shells beyond the profile contribute
exactly zero, not a small number and not a slowly growing tail. A tail summed
over every residue would pull the particle toward the chain's centroid and
fabricate the encapsulation the analysis is meant to detect, so the profile is
built to stop rather than to decay.

Sublattice parity means each residue only ever sees half the shells: residues on
the opposite sublattice to the particle sit at odd squared distances (1, 3, 5,
...) and residues sharing its sublattice at even ones (0, 2, 4, ...). A profile
covering both parities therefore describes the whole neighbourhood.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from constants import (
    DEFAULT_HARDCORE_ENERGY,
    DEFAULT_SHELL_PROFILE,
)
from logger import get_logger

if TYPE_CHECKING:
    from numpy.typing import NDArray

logger = get_logger()


class RadialPotential:
    """A finite-range, shell-resolved interaction profile.

    Attributes:
        profile (dict[float, float]): Multiplier applied to a residue's contact
            energy at each squared lattice distance. Shells absent from the
            mapping contribute exactly nothing.
        hardcore_energy (float): Residue-independent penalty for the particle
            occupying a residue's own lattice site.

    """

    def __init__(
        self,
        profile: dict[float, float] | None = None,
        hardcore_energy: float = DEFAULT_HARDCORE_ENERGY,
    ) -> None:
        """Initialise a radial potential.

        Args:
            profile (dict[float, float] | None, optional): Squared distance to
                energy multiplier. Defaults to DEFAULT_SHELL_PROFILE.
            hardcore_energy (float, optional): Penalty for site overlap. Should
                be positive and large compared with the contact energies.
                Defaults to DEFAULT_HARDCORE_ENERGY.

        Raises:
            ValueError: If *hardcore_energy* or any profile entry is not finite,
                or if a shell is negative.

        """
        if not math.isfinite(hardcore_energy):
            msg: str = f"hardcore_energy must be finite, got {hardcore_energy!r}"
            raise ValueError(msg)

        resolved: dict[float, float] = dict(
            DEFAULT_SHELL_PROFILE if profile is None else profile
        )

        for shell, weight in resolved.items():
            if shell < 0:
                msg: str = f"shell must be non-negative, got {shell!r}"
                raise ValueError(msg)
            if not math.isfinite(weight):
                msg: str = f"profile weight at shell {shell!r} must be finite"
                raise ValueError(msg)

        self.profile: dict[float, float] = resolved
        self.hardcore_energy: float = hardcore_energy

        logger.info(
            "RadialPotential created [shells=%s, hardcore=%s]",
            sorted(self.profile),
            hardcore_energy,
        )

    @classmethod
    def contact_only(
        cls, hardcore_energy: float = DEFAULT_HARDCORE_ENERGY
    ) -> RadialPotential:
        """Create a potential that acts only at touching distance.

        Reproduces the reach of the contact-flag coupling, which makes it the
        right control when comparing the two variants.

        Args:
            hardcore_energy (float, optional): Penalty for site overlap. Defaults
                to DEFAULT_HARDCORE_ENERGY.

        Returns:
            RadialPotential: A potential supported on the first shell alone.

        """
        return cls(profile={1.0: 1.0}, hardcore_energy=hardcore_energy)

    @classmethod
    def square_well(
        cls,
        depth_by_shell: dict[float, float],
        hardcore_energy: float = DEFAULT_HARDCORE_ENERGY,
    ) -> RadialPotential:
        """Create a potential from an explicit shell profile.

        Args:
            depth_by_shell (dict[float, float]): Multiplier per squared distance.
            hardcore_energy (float, optional): Penalty for site overlap. Defaults
                to DEFAULT_HARDCORE_ENERGY.

        Returns:
            RadialPotential: A potential with the requested profile.

        """
        return cls(profile=depth_by_shell, hardcore_energy=hardcore_energy)

    @classmethod
    def decaying(
        cls,
        cutoff: float = 4.0,
        falloff: float = 2.0,
        hardcore_energy: float = DEFAULT_HARDCORE_ENERGY,
    ) -> RadialPotential:
        """Create a potential that weakens with distance and stops at a cutoff.

        The weight at squared distance ``s`` is ``1 / s**(falloff / 2)``, i.e. an
        inverse power of the true separation, truncated to zero beyond *cutoff*.

        Args:
            cutoff (float, optional): Largest squared distance that still
                interacts. Defaults to 4.0.
            falloff (float, optional): Power of the distance in the denominator.
                Defaults to 2.0.
            hardcore_energy (float, optional): Penalty for site overlap. Defaults
                to DEFAULT_HARDCORE_ENERGY.

        Returns:
            RadialPotential: A truncated power-law potential.

        Raises:
            ValueError: If *cutoff* is below the first shell.

        """
        if cutoff < 1.0:
            msg: str = f"cutoff must reach at least the first shell, got {cutoff!r}"
            raise ValueError(msg)

        profile: dict[float, float] = {
            float(shell): float(shell ** (-falloff / 2.0))
            for shell in range(1, int(cutoff) + 1)
        }

        return cls(profile=profile, hardcore_energy=hardcore_energy)

    @property
    def cutoff(self) -> float:
        """float: Largest squared distance that still interacts."""
        return max(self.profile, default=0.0)

    def weight(self, squared_distance: float) -> float:
        """Return the energy multiplier at one squared distance.

        Args:
            squared_distance (float): Squared lattice distance.

        Returns:
            float: Multiplier, or 0.0 outside the profile.

        """
        for shell, value in self.profile.items():
            if math.isclose(shell, squared_distance):
                return value

        return 0.0

    def evaluate(
        self, squared_distances: NDArray[np.float64], contact_energy: float
    ) -> NDArray[np.float64]:
        """Evaluate the potential over an array of squared distances.

        Shaped for :func:`~utils.spectral_utils.build_diagonal_function`, which
        applies it to every eigenvalue of a squared-distance operator at once.

        Args:
            squared_distances (NDArray[np.float64]): Squared lattice distances.
            contact_energy (float): The residue's own interaction energy, scaled
                by the profile.

        Returns:
            NDArray[np.float64]: Energy at each distance. Exactly zero beyond the
            profile's support.

        """
        energies: NDArray[np.float64] = np.zeros_like(squared_distances)

        for shell, multiplier in self.profile.items():
            energies[np.isclose(squared_distances, shell)] = contact_energy * multiplier

        energies[np.isclose(squared_distances, 0.0)] = self.hardcore_energy

        return energies

    def __repr__(self) -> str:
        """Return a developer-readable representation of the potential."""
        shells: str = ", ".join(
            f"{shell:g}:{weight:g}" for shell, weight in sorted(self.profile.items())
        )
        return f"RadialPotential({{{shells}}}, hardcore={self.hardcore_energy:g})"
