"""End-to-end experiment runner for the ligand and external-field extensions.

Produces everything needed to judge the extension in one pass, writing to
``output/experiments/<timestamp>``:

* an apo/holo comparison with binding energy, contacts and both conformations,
  solved exactly and variationally so the two can be compared;
* an affinity sweep showing how binding responds to ligand strength;
* two field sweeps - a gradient field that reshapes the fold, and a uniform
  field that provably cannot, included as a control.

Run with::

    uv run python src/experiments.py
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

from analysis.ligand_analysis import (
    BindingResult,
    LigandAnalysis,
    SweepRecord,
    sweep_field_strength,
    sweep_ligand_affinity,
    write_records_csv,
)
from analysis.reporting import (
    plot_binding_geometry,
    plot_sweep,
    plot_vqe_convergence,
    write_binding_report,
    write_structure_xyz,
)
from constants import DEFAULT_TIMEZONE, ROOT_PROJECT_PATH
from interaction import LigandInteraction
from logger import get_logger

logger = get_logger()

DEFAULT_SEQUENCE: str = "APRLRFY"
AFFINITY_SCALES: list[float] = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0]
FIELD_STRENGTHS: list[float] = [
    -2.0,
    -1.0,
    -0.5,
    -0.25,
    -0.1,
    0.0,
    0.1,
    0.25,
    0.5,
    1.0,
    2.0,
]


def run_binding_experiment(sequence: str, output_dir: Path) -> BindingResult:
    """Compare the chain with and without a hydrophobic ligand.

    Args:
        sequence (str): Main chain sequence.
        output_dir (Path): Directory to write artefacts into.

    Returns:
        BindingResult: The comparison, already written to disk.

    """
    logger.info("Running apo/holo binding experiment for %s", sequence)

    analysis = LigandAnalysis(
        sequence=sequence,
        ligand_interaction=LigandInteraction.hp_like("H"),
        num_ligand_steps=2,
    )
    result: BindingResult = analysis.run(run_vqe=True)

    write_binding_report(result, output_dir / "binding_report.json")
    write_structure_xyz(result.apo.structure, sequence, output_dir / "apo.xyz")
    write_structure_xyz(
        result.holo.structure,
        sequence,
        output_dir / "holo.xyz",
        ligand_symbol=result.ligand.symbol,
    )
    plot_binding_geometry(result, output_dir / "binding_geometry.png")
    plot_vqe_convergence(result, output_dir / "vqe_convergence.png")

    return result


def run_affinity_sweep(sequence: str, output_dir: Path) -> list[SweepRecord]:
    """Sweep the ligand's affinity and record how binding responds.

    Args:
        sequence (str): Main chain sequence.
        output_dir (Path): Directory to write artefacts into.

    Returns:
        list[SweepRecord]: The sweep records.

    """
    logger.info("Running ligand affinity sweep for %s", sequence)

    records = sweep_ligand_affinity(sequence=sequence, energy_scales=AFFINITY_SCALES)
    write_records_csv(records, output_dir / "affinity_sweep.csv")
    plot_sweep(
        records,
        parameter_key="energy_scale",
        filepath=output_dir / "affinity_sweep.png",
        title=f"Binding vs ligand affinity ({sequence})",
        energy_key="binding_energy",
        conformation_key="holo_turns",
    )

    return records


def run_field_sweeps(
    sequence: str, output_dir: Path
) -> tuple[list[SweepRecord], list[SweepRecord]]:
    """Sweep both field modes, gradient against the uniform control.

    Args:
        sequence (str): Main chain sequence.
        output_dir (Path): Directory to write artefacts into.

    Returns:
        tuple[list[SweepRecord], list[SweepRecord]]: Gradient and
        uniform sweep records.

    """
    logger.info("Running external field sweeps for %s", sequence)

    gradient = sweep_field_strength(sequence, FIELD_STRENGTHS, gradient=True)
    uniform = sweep_field_strength(sequence, FIELD_STRENGTHS, gradient=False)

    write_records_csv(gradient, output_dir / "field_sweep_gradient.csv")
    write_records_csv(uniform, output_dir / "field_sweep_uniform.csv")

    plot_sweep(
        gradient,
        parameter_key="strength",
        filepath=output_dir / "field_sweep_gradient.png",
        title=f"Gradient field: fold responds ({sequence})",
    )
    plot_sweep(
        uniform,
        parameter_key="strength",
        filepath=output_dir / "field_sweep_uniform.png",
        title=f"Uniform field control: fold cannot respond ({sequence})",
    )

    return gradient, uniform


def summarise(
    result: BindingResult,
    affinity: list[SweepRecord],
    gradient: list[SweepRecord],
    uniform: list[SweepRecord],
    output_dir: Path,
) -> str:
    """Build a human-readable summary of every experiment.

    Args:
        result (BindingResult): Apo/holo comparison.
        affinity (list[SweepRecord]): Affinity sweep records.
        gradient (list[SweepRecord]): Gradient field sweep records.
        uniform (list[SweepRecord]): Uniform field sweep records.
        output_dir (Path): Directory the artefacts were written to.

    Returns:
        str: The summary text, also written to ``summary.txt``.

    """
    gradient_folds: set[str] = {record["turns"] for record in gradient}
    uniform_folds: set[str] = {record["turns"] for record in uniform}

    lines: list[str] = [
        f"Sequence                : {result.sequence}",
        f"Ligand                  : {result.ligand!r}",
        f"Eligible residues       : {result.ligand.eligible_bead_indices(len(result.sequence))}",
        "",
        "--- Apo vs holo (exact diagonalisation) ---",
        f"Apo  : {result.apo.num_qubits:>2} qubits, E = {result.apo.exact.energy:+.6f}, "
        f"degeneracy {result.apo.exact.degeneracy}",
        f"Holo : {result.holo.num_qubits:>2} qubits, E = {result.holo.exact.energy:+.6f}, "
        f"degeneracy {result.holo.exact.degeneracy}",
        f"Binding energy          : {result.binding_energy:+.6f}",
        f"Bound residues          : {result.bound_residues} -> {result.bound_symbols()}",
        f"Fold changed on binding : {result.conformation_changed}",
        "",
        "--- VQE vs exact ---",
        f"Apo  : VQE {result.apo.vqe_energy:+.6f}, error {result.apo.vqe_error:+.2e}, "
        f"found ground state: {result.apo.vqe_found_ground_state}",
        f"Holo : VQE {result.holo.vqe_energy:+.6f}, error {result.holo.vqe_error:+.2e}, "
        f"found ground state: {result.holo.vqe_found_ground_state}",
        "",
        "--- Ligand affinity sweep ---",
    ]

    lines.extend(
        f"  scale {record['energy_scale']:>5}  dE = {float(record['binding_energy']):+.4f}"
        f"  bound {record['bound_symbols'] or '-':<4}  fold {record['holo_turns']}"
        for record in affinity
    )

    lines.extend(
        [
            "",
            "--- External field: gradient vs uniform control ---",
            f"Gradient field produced {len(gradient_folds)} distinct fold(s) across "
            f"{len(gradient)} strengths.",
            f"Uniform field produced {len(uniform_folds)} distinct fold(s) across "
            f"{len(uniform)} strengths.",
        ]
    )

    if len(uniform_folds) == 1 and len(gradient_folds) > 1:
        lines.append(
            "Control passed: only the position-coupled field reshapes the chain."
        )
    else:
        lines.append(
            "Control INCONCLUSIVE: inspect the sweep CSVs before drawing conclusions."
        )

    for record in gradient:
        lines.append(
            f"  strength {float(record['strength']):>6.2f}  E = {float(record['energy']):+.4f}"
            f"  span = {float(record['extent_along_field']):.3f}  fold {record['turns']}"
        )

    lines.extend(["", f"Artefacts written to: {output_dir}"])

    summary: str = "\n".join(lines)
    (output_dir / "summary.txt").write_text(summary, encoding="utf-8")

    return summary


def main() -> None:
    """Run every experiment and write the artefacts plus a summary."""
    timestamp: str = datetime.now(tz=DEFAULT_TIMEZONE).strftime("%Y_%m_%d-%H_%M_%S")
    output_dir: Path = (
        ROOT_PROJECT_PATH / "output" / "experiments" / f"{timestamp}-{DEFAULT_SEQUENCE}"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    result: BindingResult = run_binding_experiment(DEFAULT_SEQUENCE, output_dir)
    affinity = run_affinity_sweep(DEFAULT_SEQUENCE, output_dir)
    gradient, uniform = run_field_sweeps(DEFAULT_SEQUENCE, output_dir)

    print(summarise(result, affinity, gradient, uniform, output_dir))


if __name__ == "__main__":
    main()
