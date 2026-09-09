import json

import numpy as np
import pytest

from analysis.ligand_analysis import (
    LigandAnalysis,
    sweep_field_strength,
    sweep_ligand_affinity,
    write_records_csv,
)
from analysis.reporting import write_binding_report, write_structure_xyz
from interaction import LigandInteraction
from particle import Ligand

SEQUENCE = "APRLRFY"


@pytest.fixture(scope="module")
def binding_result():
    analysis = LigandAnalysis(
        sequence=SEQUENCE,
        ligand_interaction=LigandInteraction.hp_like("H"),
        num_ligand_steps=2,
    )
    return analysis.run()


def test_holo_uses_more_qubits_than_apo(binding_result):
    assert binding_result.holo.num_qubits > binding_result.apo.num_qubits


def test_binding_energy_is_the_energy_difference(binding_result):
    assert np.isclose(
        binding_result.binding_energy,
        binding_result.holo.exact.energy - binding_result.apo.exact.energy,
    )


def test_hydrophobic_binding_is_favourable(binding_result):
    assert binding_result.binding_energy < 0


def test_bound_residues_are_reported_with_symbols(binding_result):
    assert binding_result.bound_residues
    assert len(binding_result.bound_symbols()) == len(binding_result.bound_residues)


def test_bound_symbols_match_the_sequence(binding_result):
    assert binding_result.bound_symbols() == [
        SEQUENCE[index] for index in binding_result.bound_residues
    ]


def test_apo_has_no_ligand_geometry(binding_result):
    assert binding_result.apo.structure.ligand_position is None


def test_holo_has_ligand_geometry(binding_result):
    assert binding_result.holo.structure.ligand_position is not None


def test_without_vqe_no_variational_result_is_recorded(binding_result):
    assert binding_result.apo.vqe_energy is None
    assert not binding_result.apo.vqe_found_ground_state
    assert binding_result.apo.vqe_error is None


def test_inert_ligand_gives_zero_binding_energy():
    analysis = LigandAnalysis(
        sequence=SEQUENCE,
        ligand_interaction=LigandInteraction.hp_like("H", energy_scale=0.0),
        num_ligand_steps=2,
    )

    assert np.isclose(analysis.run().binding_energy, 0.0)


def test_affinity_sweep_returns_one_record_per_scale():
    records = sweep_ligand_affinity(SEQUENCE, energy_scales=[0.0, 1.0])

    assert len(records) == 2
    assert [record["energy_scale"] for record in records] == [0.0, 1.0]


def test_binding_strengthens_monotonically_with_affinity():
    records = sweep_ligand_affinity(SEQUENCE, energy_scales=[0.0, 1.0, 2.0])
    energies = [float(record["binding_energy"]) for record in records]

    assert energies == sorted(energies, reverse=True)


def test_uniform_field_sweep_never_changes_the_fold():
    records = sweep_field_strength(SEQUENCE, [-2.0, 0.0, 2.0], gradient=False)

    assert len({record["turns"] for record in records}) == 1


def test_gradient_field_sweep_changes_the_fold():
    records = sweep_field_strength(SEQUENCE, [-2.0, 0.0, 2.0], gradient=True)

    assert len({record["turns"] for record in records}) > 1


def test_field_sweep_records_the_mode():
    records = sweep_field_strength(SEQUENCE, [0.0], gradient=True)

    assert records[0]["mode"] == "gradient"


def test_csv_export_writes_a_header_and_rows(tmp_path):
    records = sweep_ligand_affinity(SEQUENCE, energy_scales=[0.0, 1.0])

    path = write_records_csv(records, tmp_path / "sweep.csv")
    lines = path.read_text(encoding="utf-8").strip().splitlines()

    assert len(lines) == 3
    assert "binding_energy" in lines[0]


def test_csv_export_rejects_empty_records(tmp_path):
    with pytest.raises(ValueError, match="No records"):
        write_records_csv([], tmp_path / "empty.csv")


def test_xyz_export_includes_the_ligand(tmp_path, binding_result):
    path = write_structure_xyz(
        binding_result.holo.structure, SEQUENCE, tmp_path / "holo.xyz"
    )
    lines = path.read_text(encoding="utf-8").strip().splitlines()

    assert int(lines[0]) == len(SEQUENCE) + 1
    assert lines[-1].split()[0] == binding_result.ligand.symbol


def test_xyz_export_omits_an_absent_ligand(tmp_path, binding_result):
    path = write_structure_xyz(
        binding_result.apo.structure, SEQUENCE, tmp_path / "apo.xyz"
    )

    assert int(path.read_text(encoding="utf-8").splitlines()[0]) == len(SEQUENCE)


def test_binding_report_round_trips_through_json(tmp_path, binding_result):
    path = write_binding_report(binding_result, tmp_path / "report.json")
    payload = json.loads(path.read_text(encoding="utf-8"))

    assert payload["sequence"] == SEQUENCE
    assert np.isclose(payload["binding_energy"], binding_result.binding_energy)
    assert payload["holo"]["ligand_position"] is not None


def test_custom_ligand_can_target_a_chosen_residue():
    """F sits at index 5, so only an even-step ligand shares its sublattice."""
    analysis = LigandAnalysis(
        sequence=SEQUENCE,
        ligand_interaction=LigandInteraction.custom({"F": -3.0}, default_energy=0.0),
        ligand=Ligand(num_steps=2),
    )
    result = analysis.run()

    assert result.bound_symbols() == ["F"]


def test_a_ligand_cannot_bind_across_sublattices():
    """F is unreachable for an odd-step ligand however strong the attraction."""
    analysis = LigandAnalysis(
        sequence=SEQUENCE,
        ligand_interaction=LigandInteraction.custom({"F": -3.0}, default_energy=0.0),
        ligand=Ligand(num_steps=3),
    )
    result = analysis.run()

    assert "F" not in result.bound_symbols()
    assert np.isclose(result.binding_energy, 0.0)
