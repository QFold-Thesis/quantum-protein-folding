"""Tests for :class:`~interaction.ligand_interaction.LigandInteraction`.

Covers:
* HP_LIKE mode - via factory and direct assertions on all 20 standard residues.
* CUSTOM mode - explicit map, default fallback, partial maps.
* Factory convenience: hp_like() / custom() class-methods.
* all_energies() helper for both modes.
* Validation / error handling: bad symbols, non-finite values, wrong types.
* __repr__ sanity.
* Physical consistency checks (HH contact is the most attractive in HP-like).
"""

from __future__ import annotations

import math

import pytest

from interaction.ligand_interaction import LigandInteraction, LigandInteractionMode

# ---------------------------------------------------------------------------
# Fixtures - standard amino acid sets
# ---------------------------------------------------------------------------

# Hydrophobic residues according to the project's hp_matrix.txt
HP_HYDROPHOBIC = frozenset(["A", "C", "I", "L", "M", "F", "P", "W", "Y", "V"])
# Polar residues (remaining 10 out of 20 standard AAs)
HP_POLAR = frozenset(["R", "N", "D", "E", "Q", "G", "H", "K", "S", "T"])
ALL_STANDARD_AA = HP_HYDROPHOBIC | HP_POLAR


# ---------------------------------------------------------------------------
# Helper factories
# ---------------------------------------------------------------------------


def make_hp_h(symbol: str = "L") -> LigandInteraction:
    """Hydrophobic ligand in HP-like mode."""
    return LigandInteraction.hp_like(ligand_hp_type="H", ligand_symbol=symbol)


def make_hp_p(symbol: str = "L") -> LigandInteraction:
    """Polar ligand in HP-like mode."""
    return LigandInteraction.hp_like(ligand_hp_type="P", ligand_symbol=symbol)


def make_custom(
    energy_map: dict[str, float],
    default_energy: float = 0.0,
    symbol: str = "L",
) -> LigandInteraction:
    return LigandInteraction.custom(
        energy_map=energy_map,
        default_energy=default_energy,
        ligand_symbol=symbol,
    )


# ===========================================================================
# HP_LIKE mode
# ===========================================================================


class TestHPLikeMode:
    """Validate the HP-like interaction model for a hydrophobic/polar ligand."""

    # -- Basic attributes --

    def test_mode_is_hp_like(self):
        lig = make_hp_h()
        assert lig.mode is LigandInteractionMode.HP_LIKE

    def test_ligand_symbol_stored(self):
        lig = make_hp_h(symbol="ATP")
        assert lig.ligand_symbol == "ATP"

    def test_valid_symbols_contains_all_20_aa(self):
        lig = make_hp_h()
        assert ALL_STANDARD_AA.issubset(lig.valid_amino_acid_symbols)

    # -- Hydrophobic ligand energies --

    @pytest.mark.parametrize("aa", sorted(HP_HYDROPHOBIC))
    def test_hydrophobic_ligand_vs_hydrophobic_residue(self, aa: str):
        """H ligand + H residue must return HH energy (-1.0)."""
        lig = make_hp_h()
        assert lig.get_energy(aa) == pytest.approx(-1.0)

    @pytest.mark.parametrize("aa", sorted(HP_POLAR))
    def test_hydrophobic_ligand_vs_polar_residue(self, aa: str):
        """H ligand + P residue must return 0.0."""
        lig = make_hp_h()
        assert lig.get_energy(aa) == pytest.approx(0.0)

    # -- Polar ligand energies --

    @pytest.mark.parametrize("aa", sorted(HP_HYDROPHOBIC))
    def test_polar_ligand_vs_hydrophobic_residue(self, aa: str):
        """P ligand + H residue → 0.0."""
        lig = make_hp_p()
        assert lig.get_energy(aa) == pytest.approx(0.0)

    @pytest.mark.parametrize("aa", sorted(HP_POLAR))
    def test_polar_ligand_vs_polar_residue(self, aa: str):
        """P ligand + P residue → 0.0."""
        lig = make_hp_p()
        assert lig.get_energy(aa) == pytest.approx(0.0)

    # -- Custom energy parameters --

    def test_custom_hh_energy(self):
        lig = LigandInteraction.hp_like("H", hp_hh_energy=-3.0, hp_non_hh_energy=-0.5)
        assert lig.get_energy("A") == pytest.approx(-3.0)  # A is hydrophobic
        assert lig.get_energy("R") == pytest.approx(-0.5)  # R is polar

    def test_polar_ligand_custom_energy_never_triggers_hh(self):
        """Even with custom hp_hh_energy, a polar ligand never uses it."""
        lig = LigandInteraction.hp_like("P", hp_hh_energy=-99.0, hp_non_hh_energy=0.0)
        for aa in ALL_STANDARD_AA:
            assert lig.get_energy(aa) == pytest.approx(0.0)

    # -- Physical consistency --

    def test_hh_energy_is_more_attractive_than_non_hh(self):
        lig = make_hp_h()
        hh_energy = lig.get_energy("A")  # hydrophobic
        non_hh_energy = lig.get_energy("R")  # polar
        assert hh_energy < non_hh_energy

    # -- Error handling --

    def test_unknown_residue_raises_unsupported_symbol_error(self):
        lig = make_hp_h()
        with pytest.raises(Exception):  # noqa: B017  # UnsupportedAminoAcidSymbolError
            lig.get_energy("X")  # X not in standard HP matrix

    def test_invalid_hp_type_raises_value_error(self):
        with pytest.raises(ValueError, match="ligand_hp_type"):
            LigandInteraction.hp_like("Z")

    def test_invalid_hp_type_lowercase_raises_value_error(self):
        with pytest.raises(ValueError, match="ligand_hp_type"):
            LigandInteraction.hp_like("h")  # must be uppercase "H"

    # -- all_energies() --

    def test_all_energies_contains_all_20_residues(self):
        lig = make_hp_h()
        energies = lig.all_energies()
        assert ALL_STANDARD_AA.issubset(energies.keys())

    def test_all_energies_values_are_finite(self):
        for lig in (make_hp_h(), make_hp_p()):
            for val in lig.all_energies().values():
                assert math.isfinite(val)

    def test_all_energies_hydrophobic_ligand_sum(self):
        """HH count * (-1.0) should equal sum of negative energies."""
        lig = make_hp_h()
        energies = lig.all_energies()
        hh_count = sum(1 for aa in HP_HYDROPHOBIC if aa in energies)
        expected_sum = hh_count * (-1.0)
        actual_sum = sum(e for e in energies.values() if e < 0)
        assert actual_sum == pytest.approx(expected_sum)

    # -- repr --

    def test_repr_contains_mode(self):
        lig = make_hp_h()
        assert "HP_LIKE" in repr(lig)

    def test_repr_contains_hp_type(self):
        lig_h = make_hp_h()
        assert "'H'" in repr(lig_h)
        lig_p = make_hp_p()
        assert "'P'" in repr(lig_p)


# ===========================================================================
# CUSTOM mode
# ===========================================================================


class TestCustomMode:
    """Validate the custom energy-map interaction model."""

    # -- Basic attributes --

    def test_mode_is_custom(self):
        lig = make_custom({"A": -1.0})
        assert lig.mode is LigandInteractionMode.CUSTOM

    def test_ligand_symbol_stored(self):
        lig = make_custom({"A": -1.0}, symbol="XYZ")
        assert lig.ligand_symbol == "XYZ"

    # -- Energy lookup --

    def test_explicit_energy_returned(self):
        lig = make_custom({"A": -2.5, "K": -0.3})
        assert lig.get_energy("A") == pytest.approx(-2.5)
        assert lig.get_energy("K") == pytest.approx(-0.3)

    def test_default_energy_for_missing_residue(self):
        lig = make_custom({"A": -2.5}, default_energy=0.0)
        assert lig.get_energy("G") == pytest.approx(0.0)  # G not in map

    def test_non_zero_default_energy(self):
        lig = make_custom({"A": -1.0}, default_energy=0.5)
        assert lig.get_energy("R") == pytest.approx(0.5)

    def test_negative_default_energy(self):
        lig = make_custom({}, default_energy=-3.0)
        for aa in ALL_STANDARD_AA:
            assert lig.get_energy(aa) == pytest.approx(-3.0)

    def test_empty_energy_map_all_default(self):
        lig = make_custom({}, default_energy=-1.5)
        assert lig.get_energy("A") == pytest.approx(-1.5)
        assert lig.get_energy("K") == pytest.approx(-1.5)

    def test_energy_map_with_all_20_aa(self):
        energy_map = {aa: float(i) * -0.1 for i, aa in enumerate(ALL_STANDARD_AA)}
        lig = make_custom(energy_map)
        for aa, expected in energy_map.items():
            assert lig.get_energy(aa) == pytest.approx(expected)

    def test_positive_energy_allowed(self):
        """Custom energies are not restricted to be negative."""
        lig = make_custom({"R": +2.0, "K": +1.5})
        assert lig.get_energy("R") == pytest.approx(2.0)
        assert lig.get_energy("K") == pytest.approx(1.5)

    def test_zero_energy_allowed(self):
        lig = make_custom({"A": 0.0})
        assert lig.get_energy("A") == pytest.approx(0.0)

    # -- valid_amino_acid_symbols introspection --

    def test_valid_symbols_contains_mapped_residues(self):
        lig = make_custom({"A": -1.0, "K": -0.5})
        assert "A" in lig.valid_amino_acid_symbols
        assert "K" in lig.valid_amino_acid_symbols

    def test_valid_symbols_does_not_contain_unmapped_residue(self):
        """Unmapped residues are valid (default fallback) but not in the set."""
        lig = make_custom({"A": -1.0})
        assert "G" not in lig.valid_amino_acid_symbols  # not explicitly mapped

    # -- Validation errors --

    def test_non_finite_default_energy_raises_value_error(self):
        with pytest.raises(ValueError):
            make_custom({}, default_energy=float("inf"))

    def test_nan_default_energy_raises_value_error(self):
        with pytest.raises(ValueError):
            make_custom({}, default_energy=float("nan"))

    def test_non_finite_map_value_raises_value_error(self):
        with pytest.raises(ValueError):
            make_custom({"A": float("inf")})

    def test_nan_map_value_raises_value_error(self):
        with pytest.raises(ValueError):
            make_custom({"A": float("nan")})

    def test_non_dict_energy_map_raises_type_error(self):
        with pytest.raises(TypeError, match="energy_map"):
            LigandInteraction.custom(energy_map=[("A", -1.0)])  # type: ignore[arg-type]

    # -- all_energies() --

    def test_all_energies_returns_explicit_map(self):
        energy_map = {"A": -2.5, "K": -0.3}
        lig = make_custom(energy_map)
        assert lig.all_energies() == pytest.approx(energy_map)

    def test_all_energies_empty_map(self):
        lig = make_custom({})
        assert lig.all_energies() == {}

    def test_all_energies_does_not_include_default(self):
        """Default-energy residues should not appear in all_energies()."""
        lig = make_custom({"A": -1.0}, default_energy=-99.0)
        assert "G" not in lig.all_energies()

    # -- repr --

    def test_repr_contains_mode(self):
        lig = make_custom({"A": -1.0})
        assert "CUSTOM" in repr(lig)

    def test_repr_contains_explicit_residue_count(self):
        lig = make_custom({"A": -1.0, "K": -0.5})
        assert "2" in repr(lig)

    def test_repr_contains_default_energy(self):
        lig = make_custom({"A": -1.0}, default_energy=-3.0)
        assert "-3.0" in repr(lig)


# ===========================================================================
# Common / cross-mode tests
# ===========================================================================


class TestCommon:
    """Tests that apply to both interaction modes."""

    def test_empty_ligand_symbol_raises_value_error(self):
        with pytest.raises(ValueError, match="ligand_symbol"):
            LigandInteraction(
                mode=LigandInteractionMode.CUSTOM,
                ligand_symbol="",
                energy_map={},
            )

    def test_wrong_mode_type_raises_type_error(self):
        with pytest.raises(TypeError, match="mode"):
            LigandInteraction(
                mode="hp_like",  # type: ignore[arg-type]
                ligand_symbol="L",
            )

    def test_hp_like_and_custom_return_ligand_interaction_instances(self):
        assert isinstance(make_hp_h(), LigandInteraction)
        assert isinstance(make_custom({}), LigandInteraction)

    def test_two_independent_instances(self):
        lig_a = make_hp_h(symbol="A")
        lig_b = make_hp_p(symbol="B")
        assert lig_a.ligand_symbol != lig_b.ligand_symbol
        assert lig_a.mode is lig_b.mode  # both HP_LIKE

    def test_custom_default_ligand_symbol_is_L(self):  # noqa: N802
        lig = LigandInteraction.hp_like("H")
        assert lig.ligand_symbol == "L"

        lig2 = LigandInteraction.custom({})
        assert lig2.ligand_symbol == "L"
