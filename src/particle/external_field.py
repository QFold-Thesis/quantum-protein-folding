"""External field representation for quantum protein folding.

This module provides the :class:`ExternalField` class, which models a
position-dependent external interaction field acting on each bead of the
peptide chain on the lattice.

The field assigns an energy value to every lattice node and supports two
initialisation modes:

* **Uniform** - every node shares the same constant energy value.  Useful
  for testing whether a homogeneous background field shifts or stabilises
  particular conformations.
* **Non-uniform** - the caller supplies an explicit ``{coords: energy}``
  mapping.  Nodes absent from the mapping fall back to a configurable
  default energy (typically 0.0), so the dict only needs to list *special*
  positions.

Lattice coordinates are represented as :class:`tuple[int, ...]` of arbitrary
dimensionality (2-D square, 3-D cubic, tetrahedral …).  The class is
intentionally agnostic about the specific lattice geometry so that it can be
reused for both the simple HP model and the richer diamond/FCC lattice used by
the main quantum solver.

Example usage::

    # Uniform field with strength -0.5 at every node
    field = ExternalField.uniform(strength=-0.5)
    e = field.get_energy((1, 2))          # → -0.5

    # Non-uniform field: only nodes near the origin are special
    field = ExternalField.non_uniform(
        energy_map={(0, 0, 0): -2.0, (1, 0, 0): -1.0},
        default_energy=0.0,
    )
    e = field.get_energy((0, 0, 0))       # → -2.0
    e = field.get_energy((99, 99, 99))    # → 0.0  (default)
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

from enums import FieldMode
from logger import get_logger

if TYPE_CHECKING:
    pass

logger = get_logger()


# Type alias used throughout the module.
LatticeCoords = tuple[int, ...]


class ExternalField:
    """Position-dependent external field acting on lattice nodes.

    Stores an energy value for each lattice coordinate and exposes a single
    :meth:`get_energy` query method.  Two factory class-methods —
    :meth:`uniform` and :meth:`non_uniform` — are the intended construction
    paths; direct instantiation via ``__init__`` is also supported for
    advanced use-cases.

    Attributes:
        mode (FieldMode): Initialisation mode (UNIFORM or NON_UNIFORM).
        default_energy (float): Energy returned for coordinates that are not
            explicitly listed in the internal map.  For UNIFORM fields this
            equals the uniform strength.

    """

    def __init__(
        self,
        mode: FieldMode,
        energy_map: dict[LatticeCoords, float],
        default_energy: float,
    ) -> None:
        """Initialise an :class:`ExternalField` directly.

        Prefer the factory methods :meth:`uniform` and :meth:`non_uniform`
        over calling this constructor directly.

        Args:
            mode (FieldMode): Selects the initialisation mode that was used to
                build *energy_map*.  Stored for introspection and logging.
            energy_map (dict[LatticeCoords, float]): Explicit coordinate →
                energy mapping.  For UNIFORM fields this is typically empty
                (the default covers everything); for NON_UNIFORM fields it
                lists only the *special* positions.
            default_energy (float): Fallback energy for any coordinate absent
                from *energy_map*.

        Raises:
            TypeError: If *energy_map* is not a :class:`dict`.
            ValueError: If *default_energy* is ``NaN`` or infinite.

        """
        if not isinstance(energy_map, dict):
            msg = f"energy_map must be a dict, got {type(energy_map).__name__!r}"
            raise TypeError(msg)
        if not math.isfinite(default_energy):
            msg = f"default_energy must be a finite number, got {default_energy!r}"
            raise ValueError(msg)

        self.mode: FieldMode = mode
        self.default_energy: float = default_energy
        self._energy_map: dict[LatticeCoords, float] = dict(energy_map)

        logger.info(
            "ExternalField created [mode=%s, explicit_nodes=%d, default_energy=%s]",
            mode.name,
            len(self._energy_map),
            default_energy,
        )

    @classmethod
    def uniform(cls, strength: float) -> ExternalField:
        """Create a uniform external field with a constant energy at every node.

        Args:
            strength (float): Energy value assigned to every lattice node.

        Returns:
            ExternalField: A field instance in UNIFORM mode.

        Raises:
            ValueError: If *strength* is ``NaN`` or infinite.

        Example::

            field = ExternalField.uniform(strength=-1.0)
            assert field.get_energy((0, 0, 0)) == -1.0
            assert field.get_energy((3, 7, -2)) == -1.0

        """
        if not math.isfinite(strength):
            msg = f"strength must be a finite number, got {strength!r}"
            raise ValueError(msg)

        logger.debug("Creating uniform ExternalField with strength=%s", strength)
        return cls(
            mode=FieldMode.UNIFORM,
            energy_map={},
            default_energy=strength,
        )

    @classmethod
    def non_uniform(
        cls,
        energy_map: dict[LatticeCoords, float],
        *,
        default_energy: float = 0.0,
    ) -> ExternalField:
        """Create a non-uniform external field from an explicit coordinate map.

        Nodes absent from *energy_map* silently return *default_energy*,
        so the dict only needs to enumerate the *active* lattice positions.

        Args:
            energy_map (dict[LatticeCoords, float]): Mapping from lattice
                coordinates (arbitrary-length integer tuples) to energy values.
            default_energy (float, optional): Energy returned for coordinates
                not listed in *energy_map*.  Defaults to 0.0.

        Returns:
            ExternalField: A field instance in NON_UNIFORM mode.

        Raises:
            TypeError: If *energy_map* is not a :class:`dict`.
            ValueError: If *default_energy* or any value inside *energy_map*
                is ``NaN`` or infinite.

        Example::

            field = ExternalField.non_uniform(
                {(0, 0): -2.0, (1, 0): -1.5},
                default_energy=0.0,
            )
            assert field.get_energy((0, 0)) == -2.0
            assert field.get_energy((5, 5)) == 0.0

        """
        if not isinstance(energy_map, dict):
            msg = f"energy_map must be a dict, got {type(energy_map).__name__!r}"
            raise TypeError(msg)

        for coords, value in energy_map.items():
            if not math.isfinite(value):
                msg = (
                    f"energy_map contains a non-finite value {value!r} "
                    f"at coordinates {coords!r}"
                )
                raise ValueError(msg)

        logger.debug(
            "Creating non-uniform ExternalField with %d explicit node(s), default_energy=%s",
            len(energy_map),
            default_energy,
        )
        return cls(
            mode=FieldMode.NON_UNIFORM,
            energy_map=energy_map,
            default_energy=default_energy,
        )

    def get_energy(self, lattice_coords: LatticeCoords) -> float:
        """Return the field energy at *lattice_coords*.

        For UNIFORM fields the uniform strength is always returned (no dict
        look-up).  For NON_UNIFORM fields the explicit map is consulted first;
        if the coordinate is absent, *default_energy* is returned.

        Args:
            lattice_coords (LatticeCoords): Integer tuple identifying a node
                in the lattice (e.g. ``(x, y)`` or ``(x, y, z)``).

        Returns:
            float: Energy value at the requested coordinate.

        Raises:
            TypeError: If *lattice_coords* is not a :class:`tuple`.

        Example::

            field = ExternalField.uniform(-0.5)
            assert field.get_energy((0, 0, 0)) == -0.5

        """
        if not isinstance(lattice_coords, tuple):
            msg = (
                f"lattice_coords must be a tuple of ints, "
                f"got {type(lattice_coords).__name__!r}"
            )
            raise TypeError(msg)

        # For UNIFORM fields we still consult the explicit map first so that
        # individual nodes overwritten via set_energy() take precedence over
        # the uniform default.  If the node has not been explicitly set we
        # fall back to default_energy (which equals the uniform strength).
        energy: float = self._energy_map.get(lattice_coords, self.default_energy)
        logger.debug("ExternalField.get_energy(%s) -> %s", lattice_coords, energy)
        return energy

    def set_energy(self, lattice_coords: LatticeCoords, energy: float) -> None:
        """Explicitly set (or overwrite) the energy at *lattice_coords*.

        Works for both field modes: calling this on a UNIFORM field
        effectively converts the queried node to a *special* position
        (the mode label stays UNIFORM, but the map entry takes precedence
        for that coordinate only).

        Args:
            lattice_coords (LatticeCoords): Target node.
            energy (float): New energy value to assign.

        Raises:
            TypeError: If *lattice_coords* is not a :class:`tuple`.
            ValueError: If *energy* is ``NaN`` or infinite.

        """
        if not isinstance(lattice_coords, tuple):
            msg = (
                f"lattice_coords must be a tuple of ints, "
                f"got {type(lattice_coords).__name__!r}"
            )
            raise TypeError(msg)
        if not math.isfinite(energy):
            msg = f"energy must be a finite number, got {energy!r}"
            raise ValueError(msg)

        self._energy_map[lattice_coords] = energy
        logger.debug("ExternalField: set energy at %s -> %s", lattice_coords, energy)

    def nodes(self) -> dict[LatticeCoords, float]:
        """Return a *copy* of the explicit coordinate → energy map.

        For UNIFORM fields this will typically be an empty dict unless
        individual nodes were overwritten with :meth:`set_energy`.

        Returns:
            dict[LatticeCoords, float]: Shallow copy of the internal energy
                map.

        """
        return dict(self._energy_map)

    def __repr__(self) -> str:
        """Return a developer-readable string representation."""
        return (
            f"ExternalField(mode={self.mode.name}, "
            f"explicit_nodes={len(self._energy_map)}, "
            f"default_energy={self.default_energy})"
        )
