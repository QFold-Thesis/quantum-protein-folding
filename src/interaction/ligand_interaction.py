"""Ligand-protein interaction model for quantum protein folding.

A ligand (small molecule) interacts with the protein chain beads when it
occupies an adjacent lattice node.  This module provides
:class:`LigandInteraction`, which computes the pairwise energy between a
ligand and a single amino-acid residue.

Two independent modes are available (selected at construction time via
:class:`LigandInteractionMode`):

**HP_LIKE mode**
    The ligand is assigned a binary hydrophobic/polar character via the
    ``ligand_hp_type`` constructor argument (``"H"`` or ``"P"``).  The
    interaction energy is then determined exactly as in the standard HP model:

    * H-H contact  -> ``hp_hh_energy``   (default ``-1.0``)
    * any other pair → ``hp_non_hh_energy`` (default ``0.0``)

    This is the simplest possible extension of the HP model to ligands and
    requires no additional data: the same :data:`HP_INTERACTION_MATRIX_FILEPATH`
    resource that ``HPInteraction`` already uses is consulted to classify each
    amino acid as H or P.

**CUSTOM mode**
    The caller supplies an explicit ``energy_map: dict[str, float]`` that maps
    each amino-acid one-letter symbol to the ligand's contact energy with that
    residue.  This allows, for example, experimentally or computationally
    derived binding affinities to be plugged directly into the Hamiltonian
    later.  A configurable ``default_energy`` (default ``0.0``) is returned for
    any residue not found in the map.

Design decision - standalone class, not a subclass of :class:`Interaction`
--------------------------------------------------------------------------
:class:`Interaction` was designed for *symmetric* amino-acid-amino-acid pairs
loaded from a file.  Its interface ``get_energy(symbol_i, symbol_j)`` treats
both arguments identically.  A ligand interaction is *asymmetric*: the first
argument is always an amino-acid residue and the second is always the ligand.
Moreover, the CUSTOM mode carries no matrix file at all.

Forcing ``LigandInteraction`` into the ``Interaction`` hierarchy would require:

1. Passing a dummy file path to the superclass (misleading).
2. Accepting ``ligand_symbol`` as ``symbol_j`` and silently ignoring it.
3. No gain: no shared implementation is inherited from the ABC.

A standalone class avoids all of this while staying fully compatible with the
rest of the codebase (it exposes the same ``get_energy`` entry-point and
``valid_amino_acid_symbols`` attribute that callers need).

Example usage::

    from interaction.ligand_interaction import LigandInteraction, LigandInteractionMode

    # HP-like: ligand behaves as a hydrophobic molecule
    lig_int = LigandInteraction.hp_like(ligand_hp_type="H")
    energy = lig_int.get_energy("A")   # A is hydrophobic → -1.0
    energy = lig_int.get_energy("R")   # R is polar → 0.0

    # Custom: user-supplied energy vector
    lig_int2 = LigandInteraction.custom(
        energy_map={"A": -2.5, "K": -0.3},
        default_energy=0.0,
    )
    energy = lig_int2.get_energy("A")   # → -2.5
    energy = lig_int2.get_energy("G")   # → 0.0  (default)
"""

from __future__ import annotations

import math
from enum import Enum, auto
from pathlib import Path

import numpy as np

from constants import (
    HP_HH_CONTACT_ENERGY,
    HP_INTERACTION_MATRIX_FILEPATH,
    HP_NON_HH_CONTACT_ENERGY,
)
from exceptions import UnsupportedAminoAcidSymbolError
from logger import get_logger

logger = get_logger()

# Sentinel for "not provided" in constructor overloads
_MISSING = object()


class LigandInteractionMode(Enum):
    """Mode of the ligand-residue interaction model.

    Attributes:
        HP_LIKE: Ligand treated as hydrophobic (H) or polar (P); energies
            follow the classical HP contact rules.
        CUSTOM: Caller supplies an explicit residue → energy mapping.

    """

    HP_LIKE = auto()
    CUSTOM = auto()


class LigandInteraction:
    """Pairwise energy model between a ligand and protein-chain residues.

    Do not instantiate directly; use the factory class-methods
    :meth:`hp_like` and :meth:`custom`.

    Attributes:
        mode (LigandInteractionMode): Active interaction mode.
        ligand_symbol (str): Identifier symbol for this ligand instance.
        valid_amino_acid_symbols (frozenset[str]): Residue symbols that the
            model can handle.  For HP_LIKE this is the full set loaded from
            the HP matrix file.  For CUSTOM it is the set of keys explicitly
            present in ``energy_map``; residues *outside* this set are served
            by ``default_energy`` and are therefore also valid (the attribute
            is set to ``frozenset()`` to indicate "all residues accepted").

    """

    # ------------------------------------------------------------------
    # Internal constructor
    # ------------------------------------------------------------------

    def __init__(
        self,
        mode: LigandInteractionMode,
        ligand_symbol: str,
        *,
        # HP_LIKE fields
        ligand_hp_type: str | None = None,
        hydrophobic_symbols: frozenset[str] | None = None,
        hp_hh_energy: float = HP_HH_CONTACT_ENERGY,
        hp_non_hh_energy: float = HP_NON_HH_CONTACT_ENERGY,
        hp_valid_symbols: frozenset[str] | None = None,
        # CUSTOM fields
        energy_map: dict[str, float] | None = None,
        default_energy: float = 0.0,
    ) -> None:
        """Low-level constructor.  Prefer factory methods :meth:`hp_like` and :meth:`custom`."""
        if not ligand_symbol:
            msg = "ligand_symbol must be a non-empty string."
            raise ValueError(msg)
        if not isinstance(mode, LigandInteractionMode):
            msg = f"mode must be a LigandInteractionMode, got {type(mode).__name__!r}."
            raise TypeError(msg)

        self.mode: LigandInteractionMode = mode
        self.ligand_symbol: str = ligand_symbol

        if mode is LigandInteractionMode.HP_LIKE:
            if ligand_hp_type not in ("H", "P"):
                msg = f"ligand_hp_type must be 'H' or 'P', got {ligand_hp_type!r}."
                raise ValueError(msg)
            if hydrophobic_symbols is None:
                msg = "hydrophobic_symbols must be provided in HP_LIKE mode."
                raise ValueError(msg)
            if hp_valid_symbols is None:
                msg = "hp_valid_symbols must be provided in HP_LIKE mode."
                raise ValueError(msg)

            self._ligand_is_hydrophobic: bool = ligand_hp_type == "H"
            self._hydrophobic_symbols: frozenset[str] = hydrophobic_symbols
            self._hp_hh_energy: float = hp_hh_energy
            self._hp_non_hh_energy: float = hp_non_hh_energy
            self.valid_amino_acid_symbols: frozenset[str] = hp_valid_symbols

        else:  # CUSTOM
            if energy_map is None:
                energy_map = {}
            if not isinstance(energy_map, dict):
                msg = f"energy_map must be a dict, got {type(energy_map).__name__!r}."
                raise TypeError(msg)
            if not math.isfinite(default_energy):
                msg = f"default_energy must be finite, got {default_energy!r}."
                raise ValueError(msg)
            for aa, val in energy_map.items():
                if not math.isfinite(val):
                    msg = f"energy_map[{aa!r}] = {val!r} is not finite."
                    raise ValueError(msg)

            self._energy_map: dict[str, float] = dict(energy_map)
            self._default_energy: float = default_energy
            # In CUSTOM mode we accept any residue (unknown ones fall back to
            # default_energy), so valid_amino_acid_symbols is "open".
            # We expose the explicitly mapped symbols for introspection.
            self.valid_amino_acid_symbols: frozenset[str] = frozenset(energy_map)

        logger.info(
            "LigandInteraction created [mode=%s, ligand=%s]",
            self.mode.name,
            self.ligand_symbol,
        )

    # ------------------------------------------------------------------
    # Factory class-methods
    # ------------------------------------------------------------------

    @classmethod
    def hp_like(
        cls,
        ligand_hp_type: str,
        *,
        ligand_symbol: str = "L",
        hp_matrix_path: Path = HP_INTERACTION_MATRIX_FILEPATH,
        hp_hh_energy: float = HP_HH_CONTACT_ENERGY,
        hp_non_hh_energy: float = HP_NON_HH_CONTACT_ENERGY,
    ) -> LigandInteraction:
        """Create an HP-like ligand interaction model.

        The ligand is assigned a fixed HP character (``"H"`` or ``"P"``).
        Residue classifications are loaded from the same HP matrix file used
        by :class:`~interaction.hp_interaction.HPInteraction`.

        Args:
            ligand_hp_type (str): ``"H"`` if the ligand is hydrophobic,
                ``"P"`` if it is polar.
            ligand_symbol (str, optional): Identifier for the ligand.
                Defaults to ``"L"``.
            hp_matrix_path (Path, optional): Path to the HP matrix file.
                Defaults to the project-level ``HP_INTERACTION_MATRIX_FILEPATH``.
            hp_hh_energy (float, optional): Energy for a hydrophobic-hydrophobic
                contact.  Defaults to :data:`~constants.HP_HH_CONTACT_ENERGY`
                (``-1.0``).
            hp_non_hh_energy (float, optional): Energy for any non-HH contact.
                Defaults to :data:`~constants.HP_NON_HH_CONTACT_ENERGY`
                (``0.0``).

        Returns:
            LigandInteraction: Instance in :attr:`LigandInteractionMode.HP_LIKE`
                mode.

        Raises:
            ValueError: If ``ligand_hp_type`` is not ``"H"`` or ``"P"``.

        Example::

            lig = LigandInteraction.hp_like("H")
            assert lig.get_energy("A") == -1.0  # A is hydrophobic
            assert lig.get_energy("R") == 0.0   # R is polar

        """
        if ligand_hp_type not in ("H", "P"):
            msg = f"ligand_hp_type must be 'H' or 'P', got {ligand_hp_type!r}."
            raise ValueError(msg)

        hydrophobic, polar = cls._load_hp_symbols(hp_matrix_path)
        hp_valid = frozenset(hydrophobic) | frozenset(polar)

        return cls(
            mode=LigandInteractionMode.HP_LIKE,
            ligand_symbol=ligand_symbol,
            ligand_hp_type=ligand_hp_type,
            hydrophobic_symbols=frozenset(hydrophobic),
            hp_hh_energy=hp_hh_energy,
            hp_non_hh_energy=hp_non_hh_energy,
            hp_valid_symbols=hp_valid,
        )

    @classmethod
    def custom(
        cls,
        energy_map: dict[str, float],
        *,
        ligand_symbol: str = "L",
        default_energy: float = 0.0,
    ) -> LigandInteraction:
        """Create a custom ligand interaction model from an explicit energy map.

        Args:
            energy_map (dict[str, float]): Mapping from one-letter amino-acid
                symbol to interaction energy with the ligand.  Residues absent
                from the map receive ``default_energy``.
            ligand_symbol (str, optional): Identifier for the ligand.
                Defaults to ``"L"``.
            default_energy (float, optional): Fallback energy for residues not
                listed in ``energy_map``.  Must be finite.  Defaults to ``0.0``.

        Returns:
            LigandInteraction: Instance in :attr:`LigandInteractionMode.CUSTOM`
                mode.

        Raises:
            TypeError: If ``energy_map`` is not a :class:`dict`.
            ValueError: If ``default_energy`` or any value in ``energy_map``
                is non-finite.

        Example::

            lig = LigandInteraction.custom(
                energy_map={"A": -2.5, "K": -0.3},
                default_energy=0.0,
            )
            assert lig.get_energy("A") == -2.5
            assert lig.get_energy("G") == 0.0   # default

        """
        return cls(
            mode=LigandInteractionMode.CUSTOM,
            ligand_symbol=ligand_symbol,
            energy_map=energy_map,
            default_energy=default_energy,
        )

    # ------------------------------------------------------------------
    # Core query
    # ------------------------------------------------------------------

    def get_energy(self, amino_acid_symbol: str) -> float:
        """Return the interaction energy between the ligand and *amino_acid_symbol*.

        Args:
            amino_acid_symbol (str): One-letter amino-acid residue symbol.

        Returns:
            float: Interaction energy.

        Raises:
            UnsupportedAminoAcidSymbolError: In HP_LIKE mode, if
                *amino_acid_symbol* is not present in the loaded HP matrix.

        """
        if self.mode is LigandInteractionMode.HP_LIKE:
            return self._get_energy_hp_like(amino_acid_symbol)
        return self._get_energy_custom(amino_acid_symbol)

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _get_energy_hp_like(self, amino_acid_symbol: str) -> float:
        """Compute HP-like interaction energy."""
        if amino_acid_symbol not in self.valid_amino_acid_symbols:
            msg = (
                f"Amino acid symbol {amino_acid_symbol!r} not found in the HP "
                f"matrix used by LigandInteraction (mode=HP_LIKE)."
            )
            logger.error(msg)
            raise UnsupportedAminoAcidSymbolError(msg)

        residue_is_hydrophobic = amino_acid_symbol in self._hydrophobic_symbols
        if self._ligand_is_hydrophobic and residue_is_hydrophobic:
            energy = self._hp_hh_energy
        else:
            energy = self._hp_non_hh_energy

        logger.debug(
            "LigandInteraction(HP_LIKE): get_energy(%s) -> %s",
            amino_acid_symbol,
            energy,
        )
        return energy

    def _get_energy_custom(self, amino_acid_symbol: str) -> float:
        """Retrieve custom interaction energy (with default fallback)."""
        energy = self._energy_map.get(amino_acid_symbol, self._default_energy)
        logger.debug(
            "LigandInteraction(CUSTOM): get_energy(%s) -> %s", amino_acid_symbol, energy
        )
        return energy

    @staticmethod
    def _load_hp_symbols(hp_filepath: Path) -> tuple[list[str], list[str]]:
        """Load hydrophobic and polar symbol lists from an HP matrix file.

        Args:
            hp_filepath (Path): Path to the HP matrix text file.

        Returns:
            tuple[list[str], list[str]]: ``(hydrophobic, polar)`` lists.

        Raises:
            Exception: If the file cannot be read or parsed.

        """
        try:
            hp_matrix = np.loadtxt(hp_filepath, dtype=str)
            hydrophobic: list[str] = []
            polar: list[str] = []
            for line in hp_matrix:
                if line[1] == "1":
                    hydrophobic.append(line[0])
                else:
                    polar.append(line[0])
        except Exception:
            logger.exception("Error loading HP matrix for LigandInteraction")
            raise
        else:
            logger.debug(
                "LigandInteraction: loaded %d hydrophobic, %d polar symbols from %s",
                len(hydrophobic),
                len(polar),
                hp_filepath,
            )
            return hydrophobic, polar

    # ------------------------------------------------------------------
    # Convenience helpers
    # ------------------------------------------------------------------

    def all_energies(self) -> dict[str, float]:
        """Return a mapping of all *explicitly known* residue energies.

        For HP_LIKE mode this covers all residues in the HP matrix.
        For CUSTOM mode this covers only the residues listed in ``energy_map``;
        residues that fall back to ``default_energy`` are not included.

        Returns:
            dict[str, float]: ``{amino_acid_symbol: energy}`` for known residues.

        """
        if self.mode is LigandInteractionMode.HP_LIKE:
            return {
                aa: self._get_energy_hp_like(aa) for aa in self.valid_amino_acid_symbols
            }
        return dict(self._energy_map)

    def __repr__(self) -> str:
        """Return a developer-readable string representation."""
        if self.mode is LigandInteractionMode.HP_LIKE:
            hp_char = "H" if self._ligand_is_hydrophobic else "P"
            return (
                f"LigandInteraction(mode=HP_LIKE, ligand={self.ligand_symbol!r}, "
                f"hp_type={hp_char!r})"
            )
        return (
            f"LigandInteraction(mode=CUSTOM, ligand={self.ligand_symbol!r}, "
            f"explicit_residues={len(self._energy_map)}, "
            f"default_energy={self._default_energy})"
        )
