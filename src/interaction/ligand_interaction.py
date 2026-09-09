"""Ligand-residue interaction energies.

Unlike :class:`~interaction.interaction.Interaction`, which scores a *pair* of
residues, a ligand interaction scores a single residue against the one fixed
chemical identity of the ligand. Two ways of supplying those energies are
provided:

* **HP-like** - the ligand is declared to be hydrophobic or polar and the
  existing HP matrix decides each residue's energy. A hydrophobic ligand is
  attracted to hydrophobic residues and indifferent to polar ones, which makes
  it settle into the hydrophobic core the folding term is already building.

  Note that the ligand's H/P label is resolved separately from the residue
  matrix. The one-letter residue codes collide with the HP labels but carry the
  opposite meaning - histidine ("H") is polar and proline ("P") is hydrophobic -
  so feeding the ligand's label through the pairwise residue lookup would
  silently invert it.
* **Custom** - an explicit ``{residue: energy}`` map for modelling a specific
  binding preference, with a default energy for residues left unlisted.
"""

from __future__ import annotations

import math

from constants import (
    DEFAULT_LIGAND_SYMBOL,
    HP_HH_CONTACT_ENERGY,
    HP_NON_HH_CONTACT_ENERGY,
)
from enums import LigandInteractionMode
from exceptions import UnsupportedAminoAcidSymbolError
from interaction.hp_interaction import HPInteraction
from logger import get_logger

logger = get_logger()


class LigandInteraction:
    """Scores the contact energy between the ligand and a single residue.

    Attributes:
        mode (LigandInteractionMode): How energies are resolved.
        symbol (str): Symbol identifying the ligand.
        valid_symbols (set[str]): Residue symbols this model can score.

    """

    def __init__(
        self,
        mode: LigandInteractionMode,
        symbol: str = DEFAULT_LIGAND_SYMBOL,
        energy_map: dict[str, float] | None = None,
        default_energy: float = 0.0,
        hp_type: str | None = None,
        energy_scale: float = 1.0,
    ) -> None:
        """Initialise a ligand interaction model.

        Prefer the :meth:`hp_like` and :meth:`custom` factory methods.

        Args:
            mode (LigandInteractionMode): How energies are resolved.
            symbol (str, optional): Ligand symbol. Defaults to
                DEFAULT_LIGAND_SYMBOL.
            energy_map (dict[str, float] | None, optional): Explicit residue to
                energy mapping, required in CUSTOM mode. Defaults to None.
            default_energy (float, optional): Energy for residues absent from
                *energy_map*. Defaults to 0.0.
            hp_type (str | None, optional): Ligand's HP character, required in
                HP_LIKE mode. Defaults to None.
            energy_scale (float, optional): Affinity multiplier applied to every
                resolved energy. Sweeping it traces out binding regimes without
                changing the chemistry. Defaults to 1.0.

        Raises:
            ValueError: If *default_energy*, *energy_scale* or any mapped energy
                is not finite, or if the arguments required by *mode* are missing
                or invalid.

        """
        if not math.isfinite(default_energy):
            msg: str = f"default_energy must be finite, got {default_energy!r}"
            raise ValueError(msg)

        if not math.isfinite(energy_scale):
            msg: str = f"energy_scale must be finite, got {energy_scale!r}"
            raise ValueError(msg)

        self.mode: LigandInteractionMode = mode
        self.symbol: str = symbol
        self.default_energy: float = default_energy
        self.energy_scale: float = energy_scale
        self._energy_map: dict[str, float] = dict(energy_map or {})
        self._hp_interaction: HPInteraction | None = None
        self._hp_type: str | None = hp_type

        for residue, energy in self._energy_map.items():
            if not math.isfinite(energy):
                msg: str = (
                    f"energy_map holds a non-finite energy {energy!r} for residue "
                    f"{residue!r}"
                )
                raise ValueError(msg)

        if mode == LigandInteractionMode.HP_LIKE:
            self._init_hp_like(hp_type)
        else:
            self.valid_symbols: set[str] = set(self._energy_map)

        logger.info(
            "LigandInteraction %s initialised [mode=%s, residues=%d, default=%s]",
            self.symbol,
            self.mode.value,
            len(self.valid_symbols),
            self.default_energy,
        )

    def _init_hp_like(self, hp_type: str | None) -> None:
        """Prepare HP-backed energy resolution.

        Args:
            hp_type (str | None): The ligand's HP character, "H" or "P".

        Raises:
            ValueError: If *hp_type* is missing or is neither "H" nor "P".

        """
        if hp_type is None:
            msg: str = "hp_type is required when mode is HP_LIKE"
            raise ValueError(msg)

        normalised: str = hp_type.upper()
        if normalised not in {"H", "P"}:
            msg: str = f"hp_type must be either 'H' or 'P', got {hp_type!r}"
            raise ValueError(msg)

        self._hp_type = normalised
        self._hp_interaction = HPInteraction()
        self.valid_symbols = set(self._hp_interaction.valid_symbols)

    @classmethod
    def hp_like(
        cls,
        hp_type: str,
        symbol: str = DEFAULT_LIGAND_SYMBOL,
        energy_scale: float = 1.0,
    ) -> LigandInteraction:
        """Create a ligand scored through the hydrophobic/polar matrix.

        Args:
            hp_type (str): The ligand's own HP character, "H" or "P".
            symbol (str, optional): Ligand symbol. Defaults to
                DEFAULT_LIGAND_SYMBOL.
            energy_scale (float, optional): Affinity multiplier. Defaults to 1.0.

        Returns:
            LigandInteraction: A model in HP_LIKE mode.

        """
        return cls(
            mode=LigandInteractionMode.HP_LIKE,
            symbol=symbol,
            hp_type=hp_type,
            energy_scale=energy_scale,
        )

    @classmethod
    def custom(
        cls,
        energy_map: dict[str, float],
        default_energy: float = 0.0,
        symbol: str = DEFAULT_LIGAND_SYMBOL,
        energy_scale: float = 1.0,
    ) -> LigandInteraction:
        """Create a ligand scored through an explicit per-residue energy map.

        Args:
            energy_map (dict[str, float]): Residue symbol to contact energy.
            default_energy (float, optional): Energy for residues absent from
                the map. Defaults to 0.0.
            symbol (str, optional): Ligand symbol. Defaults to
                DEFAULT_LIGAND_SYMBOL.
            energy_scale (float, optional): Affinity multiplier. Defaults to 1.0.

        Returns:
            LigandInteraction: A model in CUSTOM mode.

        """
        return cls(
            mode=LigandInteractionMode.CUSTOM,
            symbol=symbol,
            energy_map=energy_map,
            default_energy=default_energy,
            energy_scale=energy_scale,
        )

    def get_energy(self, residue_symbol: str) -> float:
        """Return the contact energy between the ligand and one residue.

        Args:
            residue_symbol (str): Single-letter residue symbol.

        Returns:
            float: Contact energy. Negative values are attractive.

        Raises:
            UnsupportedAminoAcidSymbolError: If HP_LIKE mode is active and the
                residue is absent from the HP matrix.

        """
        if self.mode == LigandInteractionMode.CUSTOM:
            return self.energy_scale * self._energy_map.get(
                residue_symbol, self.default_energy
            )

        if self._hp_interaction is None or self._hp_type is None:
            msg: str = "HP_LIKE ligand interaction was not initialised correctly"
            raise UnsupportedAminoAcidSymbolError(msg)

        residue_is_hydrophobic: bool = self._hp_interaction.is_hydrophobic(
            residue_symbol
        )

        return self.energy_scale * (
            HP_HH_CONTACT_ENERGY
            if (self._hp_type == "H" and residue_is_hydrophobic)
            else HP_NON_HH_CONTACT_ENERGY
        )

    def __repr__(self) -> str:
        """Return a developer-readable representation of the interaction."""
        detail: str = (
            f"hp_type={self._hp_type!r}"
            if self.mode == LigandInteractionMode.HP_LIKE
            else f"residues={len(self._energy_map)}"
        )
        return (
            f"LigandInteraction(symbol={self.symbol!r}, "
            f"mode={self.mode.value!r}, {detail})"
        )
