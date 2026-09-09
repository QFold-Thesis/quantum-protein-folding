import pytest

from constants import HP_HH_CONTACT_ENERGY, HP_NON_HH_CONTACT_ENERGY
from enums import LigandInteractionMode
from exceptions import UnsupportedAminoAcidSymbolError
from interaction import HPInteraction, LigandInteraction


@pytest.fixture(scope="module")
def hp_matrix():
    return HPInteraction()


def test_hydrophobic_ligand_attracts_hydrophobic_residues(hp_matrix):
    ligand = LigandInteraction.hp_like("H")

    for symbol in hp_matrix.valid_symbols:
        expected = (
            HP_HH_CONTACT_ENERGY
            if hp_matrix.is_hydrophobic(symbol)
            else HP_NON_HH_CONTACT_ENERGY
        )
        assert ligand.get_energy(symbol) == expected


def test_polar_ligand_is_indifferent_to_every_residue(hp_matrix):
    ligand = LigandInteraction.hp_like("P")

    assert all(
        ligand.get_energy(symbol) == HP_NON_HH_CONTACT_ENERGY
        for symbol in hp_matrix.valid_symbols
    )


def test_ligand_label_is_not_read_as_an_amino_acid(hp_matrix):
    """Histidine is polar and proline hydrophobic, opposite to the HP labels."""
    assert not hp_matrix.is_hydrophobic("H")
    assert hp_matrix.is_hydrophobic("P")

    hydrophobic_ligand = LigandInteraction.hp_like("H")

    assert hydrophobic_ligand.get_energy("P") == HP_HH_CONTACT_ENERGY
    assert hydrophobic_ligand.get_energy("H") == HP_NON_HH_CONTACT_ENERGY


def test_custom_map_returns_listed_energies():
    ligand = LigandInteraction.custom({"A": -2.5, "K": -0.3}, default_energy=0.1)

    assert ligand.get_energy("A") == pytest.approx(-2.5)
    assert ligand.get_energy("K") == pytest.approx(-0.3)


def test_custom_map_falls_back_to_the_default():
    ligand = LigandInteraction.custom({"A": -2.5}, default_energy=0.7)

    assert ligand.get_energy("W") == pytest.approx(0.7)


def test_energy_scale_multiplies_hp_energies():
    ligand = LigandInteraction.hp_like("H", energy_scale=3.0)

    assert ligand.get_energy("P") == pytest.approx(3.0 * HP_HH_CONTACT_ENERGY)


def test_energy_scale_multiplies_custom_energies():
    ligand = LigandInteraction.custom({"A": -2.0}, energy_scale=0.5)

    assert ligand.get_energy("A") == pytest.approx(-1.0)


def test_zero_scale_makes_the_ligand_inert():
    ligand = LigandInteraction.hp_like("H", energy_scale=0.0)

    assert ligand.get_energy("P") == pytest.approx(0.0)


def test_factories_set_the_mode():
    assert LigandInteraction.hp_like("H").mode is LigandInteractionMode.HP_LIKE
    assert LigandInteraction.custom({}).mode is LigandInteractionMode.CUSTOM


def test_hp_type_is_case_insensitive():
    assert LigandInteraction.hp_like("h").get_energy("P") == HP_HH_CONTACT_ENERGY


@pytest.mark.parametrize("hp_type", ["X", "HP", ""])
def test_invalid_hp_type_is_rejected(hp_type):
    with pytest.raises(ValueError, match="hp_type"):
        LigandInteraction.hp_like(hp_type)


def test_missing_hp_type_is_rejected():
    with pytest.raises(ValueError, match="hp_type is required"):
        LigandInteraction(mode=LigandInteractionMode.HP_LIKE)


@pytest.mark.parametrize("energy", [float("nan"), float("inf")])
def test_non_finite_energies_are_rejected(energy):
    with pytest.raises(ValueError, match="finite"):
        LigandInteraction.custom({"A": energy})


def test_non_finite_scale_is_rejected():
    with pytest.raises(ValueError, match="energy_scale"):
        LigandInteraction.hp_like("H", energy_scale=float("nan"))


def test_unknown_residue_is_rejected_in_hp_mode():
    with pytest.raises(UnsupportedAminoAcidSymbolError):
        LigandInteraction.hp_like("H").get_energy("Z")
