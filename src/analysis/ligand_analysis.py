"""Protein-ligand binding analysis.

Runs the same chain twice - once on its own (*apo*) and once sharing the lattice
with a ligand (*holo*) - and reports what changed. Because the ligand is coupled
through the actual lattice geometry rather than through a constant offset, the
comparison is meaningful: the fold can rearrange to accommodate the ligand, and
the binding energy measures a real preference rather than a bookkeeping shift.

Each Hamiltonian is solved twice as well: exactly, by enumerating the diagonal,
and variationally with the project's VQE setup. The exact answer is the yardstick
that says whether the variational run actually found the ground state.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np

from analysis.exact_solver import (
    ExactSolution,
    compress_with_layout,
    expand_bitstring,
    solve_exactly,
)
from analysis.particle_dynamics import ParticleEnvironment, analyse_environment
from analysis.structure_decoder import DecodedStructure, decode_structure
from builder import HamiltonianBuilder
from constants import DEFAULT_LIGAND_STEPS, EMPTY_SIDECHAIN_PLACEHOLDER
from enums import LigandCoupling
from interaction import LigandInteraction
from logger import get_logger
from particle import ExternalField, Ligand, RadialPotential
from utils.setup_utils import (
    run_vqe_optimization,
    setup_folding_system,
    setup_vqe_optimization,
)

if TYPE_CHECKING:
    from pathlib import Path

    from numpy.typing import NDArray
    from qiskit.quantum_info import SparsePauliOp

logger = get_logger()

SweepRecord = dict[str, Any]
"""One row of a parameter sweep, keyed by column name."""


@dataclass(frozen=True)
class SystemOutcome:
    """Result of solving one Hamiltonian both exactly and variationally.

    Attributes:
        label (str): Human-readable name of the system, e.g. "apo" or "holo".
        num_qubits (int): Width of the compressed register that was solved.
        exact (ExactSolution): Exact ground state of the compressed Hamiltonian.
        structure (DecodedStructure): Geometry of the exact ground state.
        vqe_energy (float | None): Best energy the VQE reached, if it was run.
        vqe_bitstring (str | None): Bitstring behind *vqe_energy*.
        vqe_structure (DecodedStructure | None): Geometry the VQE settled on.
        vqe_energies (list[float]): Energy at each VQE evaluation.

    """

    label: str
    num_qubits: int
    exact: ExactSolution
    structure: DecodedStructure
    vqe_energy: float | None = None
    vqe_bitstring: str | None = None
    vqe_structure: DecodedStructure | None = None
    vqe_energies: list[float] = field(default_factory=list)

    @property
    def vqe_found_ground_state(self) -> bool:
        """bool: Whether the VQE landed on the exact ground state energy."""
        if self.vqe_energy is None:
            return False
        return bool(np.isclose(self.vqe_energy, self.exact.energy, atol=1e-6))

    @property
    def vqe_error(self) -> float | None:
        """Float | None: Gap between the VQE energy and the exact ground state."""
        if self.vqe_energy is None:
            return None
        return float(self.vqe_energy - self.exact.energy)


@dataclass(frozen=True)
class BindingResult:
    """Comparison of a chain with and without its ligand.

    Attributes:
        sequence (str): The main chain sequence studied.
        apo (SystemOutcome): The chain solved on its own.
        holo (SystemOutcome): The chain solved with the ligand present.
        ligand (Ligand): The ligand that was placed on the lattice.
        bound_residues (list[int]): Residues the ligand actually touches in the
            holo ground state.
        conformation_changed (bool): Whether binding altered the backbone turns.

    """

    sequence: str
    apo: SystemOutcome
    holo: SystemOutcome
    ligand: Ligand
    bound_residues: list[int]
    conformation_changed: bool

    @property
    def binding_energy(self) -> float:
        """float: Energy released on binding, negative when binding is favourable.

        Defined as ``E_holo - E_apo``, so it folds together the ligand contact
        energy and any strain the chain takes on to accommodate it.
        """
        return float(self.holo.exact.energy - self.apo.exact.energy)

    def bound_symbols(self) -> list[str]:
        """Return the residue letters the ligand binds.

        Returns:
            list[str]: One-letter symbols of the contacted residues.

        """
        return [self.sequence[index] for index in self.bound_residues]


class LigandAnalysis:
    """Builds, solves and compares the apo and holo forms of one chain.

    Attributes:
        sequence (str): Main chain sequence under study.
        ligand (Ligand): Ligand placed on the lattice.
        ligand_interaction (LigandInteraction): Ligand-residue energies.
        external_field (ExternalField | None): Optional field applied to both
            forms so the comparison stays fair.

    """

    def __init__(
        self,
        sequence: str,
        ligand_interaction: LigandInteraction | None = None,
        ligand: Ligand | None = None,
        num_ligand_steps: int = DEFAULT_LIGAND_STEPS,
        external_field: ExternalField | None = None,
        coupling: LigandCoupling = LigandCoupling.CONTACT_FLAG,
        radial_potential: RadialPotential | None = None,
    ) -> None:
        """Initialise the analysis.

        Args:
            sequence (str): Main chain sequence.
            ligand_interaction (LigandInteraction | None, optional): Ligand
                energies. Defaults to a hydrophobic HP-like ligand.
            ligand (Ligand | None, optional): Ligand to place. Defaults to one
                built from *num_ligand_steps*.
            num_ligand_steps (int, optional): Lattice steps encoding the ligand
                position, used when *ligand* is not given. Defaults to
                DEFAULT_LIGAND_STEPS.
            external_field (ExternalField | None, optional): Field applied to
                both forms. Defaults to None.
            coupling (LigandCoupling, optional): How the particle couples to the
                chain. Defaults to LigandCoupling.CONTACT_FLAG.
            radial_potential (RadialPotential | None, optional): Shell profile
                for RADIAL_POTENTIAL coupling. Defaults to the standard well.

        """
        self.sequence: str = sequence
        self.ligand: Ligand = ligand or Ligand(num_steps=num_ligand_steps)
        self.ligand_interaction: LigandInteraction = (
            ligand_interaction or LigandInteraction.hp_like("H")
        )
        self.external_field: ExternalField | None = external_field
        self.coupling: LigandCoupling = coupling
        self.radial_potential: RadialPotential = radial_potential or RadialPotential()

        side_chain: str = EMPTY_SIDECHAIN_PLACEHOLDER * len(sequence)
        (
            self._protein,
            self._interaction,
            self._contact_map,
            self._distance_map,
        ) = setup_folding_system(main_chain=sequence, side_chain=side_chain)

    def _build(self, *, with_ligand: bool) -> HamiltonianBuilder:
        """Create a builder for one of the two forms.

        Args:
            with_ligand (bool): Whether to include the ligand registers.

        Returns:
            HamiltonianBuilder: Configured builder.

        """
        return HamiltonianBuilder(
            protein=self._protein,
            interaction=self._interaction,
            distance_map=self._distance_map,
            contact_map=self._contact_map,
            external_field=self.external_field,
            ligand=self.ligand if with_ligand else None,
            ligand_interaction=self.ligand_interaction if with_ligand else None,
            coupling=self.coupling,
            radial_potential=self.radial_potential,
        )

    def solve(
        self, label: str, *, with_ligand: bool, run_vqe: bool = False
    ) -> SystemOutcome:
        """Build and solve one form of the system.

        Args:
            label (str): Name to record the outcome under.
            with_ligand (bool): Whether the ligand is present.
            run_vqe (bool, optional): Also run the variational solver and compare
                it against the exact answer. Defaults to False.

        Returns:
            SystemOutcome: Exact solution, decoded geometry and any VQE result.

        """
        builder: HamiltonianBuilder = self._build(with_ligand=with_ligand)
        hamiltonian: SparsePauliOp = builder.sum_hamiltonians()
        compressed, kept = compress_with_layout(hamiltonian)

        logger.info(
            "Solving %s system: %d qubits compressed to %d",
            label,
            hamiltonian.num_qubits,
            compressed.num_qubits,
        )

        exact: ExactSolution = solve_exactly(compressed)
        structure: DecodedStructure = self._decode(
            exact.bitstring, kept, hamiltonian, builder, with_ligand=with_ligand
        )

        if not run_vqe:
            return SystemOutcome(
                label=label,
                num_qubits=int(compressed.num_qubits),
                exact=exact,
                structure=structure,
            )

        vqe, _counts, values = setup_vqe_optimization(
            num_qubits=int(compressed.num_qubits)
        )
        raw = run_vqe_optimization(vqe=vqe, hamiltonian=compressed)

        best = raw.best_measurement or {}
        vqe_bitstring: str = str(best.get("bitstring", ""))
        vqe_energy: float = float(np.real(best.get("value", np.nan)))
        vqe_structure: DecodedStructure = self._decode(
            vqe_bitstring, kept, hamiltonian, builder, with_ligand=with_ligand
        )

        return SystemOutcome(
            label=label,
            num_qubits=int(compressed.num_qubits),
            exact=exact,
            structure=structure,
            vqe_energy=vqe_energy,
            vqe_bitstring=vqe_bitstring,
            vqe_structure=vqe_structure,
            vqe_energies=list(values),
        )

    def _decode(
        self,
        bitstring: str,
        kept_qubits: list[int],
        hamiltonian: SparsePauliOp,
        builder: HamiltonianBuilder,
        *,
        with_ligand: bool,
    ) -> DecodedStructure:
        """Expand a compressed bitstring and read the geometry out of it.

        Args:
            bitstring (str): State measured on the compressed register.
            kept_qubits (list[int]): Qubits that survived compression.
            hamiltonian (SparsePauliOp): The uncompressed Hamiltonian.
            builder (HamiltonianBuilder): Builder holding the register layout.
            with_ligand (bool): Whether ligand registers are present.

        Returns:
            DecodedStructure: The decoded conformation.

        """
        full: str = expand_bitstring(
            bitstring, kept_qubits, int(hamiltonian.num_qubits)
        )

        return decode_structure(
            bitstring=full,
            chain_length=len(self.sequence),
            ligand=self.ligand if with_ligand else None,
            ligand_walk_offset=builder.ligand_walk_offset,
            ligand_contact_offset=builder.contact_flag_offset,
        )

    def run(self, *, run_vqe: bool = False) -> BindingResult:
        """Solve both forms and compare them.

        Args:
            run_vqe (bool, optional): Also run the variational solver on each
                form. Defaults to False.

        Returns:
            BindingResult: The apo and holo outcomes plus the binding summary.

        """
        apo: SystemOutcome = self.solve("apo", with_ligand=False, run_vqe=run_vqe)
        holo: SystemOutcome = self.solve("holo", with_ligand=True, run_vqe=run_vqe)

        apo_turns: list[int] = [turn.value for turn in apo.structure.turns]
        holo_turns: list[int] = [turn.value for turn in holo.structure.turns]

        result = BindingResult(
            sequence=self.sequence,
            apo=apo,
            holo=holo,
            ligand=self.ligand,
            bound_residues=(
                holo.structure.realised_contacts
                if self.coupling == LigandCoupling.CONTACT_FLAG
                else holo.structure.lattice_neighbours
            ),
            conformation_changed=apo_turns != holo_turns,
        )

        logger.info(
            "Binding analysis for %s: dE=%.6f, bound residues=%s, fold changed=%s",
            self.sequence,
            result.binding_energy,
            result.bound_residues,
            result.conformation_changed,
        )
        return result


def sweep_ligand_affinity(
    sequence: str,
    energy_scales: list[float],
    hp_type: str = "H",
    num_ligand_steps: int = DEFAULT_LIGAND_STEPS,
) -> list[SweepRecord]:
    """Trace binding behaviour as the ligand's affinity is varied.

    Args:
        sequence (str): Main chain sequence.
        energy_scales (list[float]): Affinity multipliers to evaluate.
        hp_type (str, optional): The ligand's HP character. Defaults to "H".
        num_ligand_steps (int, optional): Lattice steps encoding the ligand's
            position. Defaults to DEFAULT_LIGAND_STEPS.

    Returns:
        list[SweepRecord]: One record per affinity, ready for CSV export.

    """
    records: list[SweepRecord] = []

    for scale in energy_scales:
        analysis = LigandAnalysis(
            sequence=sequence,
            ligand_interaction=LigandInteraction.hp_like(
                hp_type=hp_type, energy_scale=scale
            ),
            num_ligand_steps=num_ligand_steps,
        )
        result: BindingResult = analysis.run()

        records.append(
            {
                "energy_scale": scale,
                "apo_energy": result.apo.exact.energy,
                "holo_energy": result.holo.exact.energy,
                "binding_energy": result.binding_energy,
                "bound_residues": " ".join(map(str, result.bound_residues)),
                "bound_symbols": " ".join(result.bound_symbols()),
                "conformation_changed": result.conformation_changed,
                "apo_turns": " ".join(
                    str(turn.value) for turn in result.apo.structure.turns
                ),
                "holo_turns": " ".join(
                    str(turn.value) for turn in result.holo.structure.turns
                ),
                "holo_degeneracy": result.holo.exact.degeneracy,
            }
        )
        logger.info(
            "Affinity %.3f -> dE=%.4f, bound=%s",
            scale,
            result.binding_energy,
            result.bound_residues,
        )

    return records


def sweep_field_strength(
    sequence: str,
    strengths: list[float],
    field_factory: type[ExternalField] = ExternalField,
    *,
    gradient: bool = True,
) -> list[SweepRecord]:
    """Trace how an external field reshapes the fold.

    Running this with ``gradient=False`` is the control experiment: a uniform
    field shifts every energy identically, so the fold must stay put no matter
    how strong it gets.

    Args:
        sequence (str): Main chain sequence.
        strengths (list[float]): Field magnitudes to evaluate.
        field_factory (type[ExternalField], optional): Field class to build with.
            Defaults to ExternalField.
        gradient (bool, optional): Whether to couple the field to bead positions.
            Defaults to True.

    Returns:
        list[SweepRecord]: One record per field strength.

    """
    records: list[SweepRecord] = []
    side_chain: str = EMPTY_SIDECHAIN_PLACEHOLDER * len(sequence)
    protein, interaction, contact_map, distance_map = setup_folding_system(
        main_chain=sequence, side_chain=side_chain
    )

    for strength in strengths:
        field: ExternalField = (
            field_factory.gradient(strength)
            if gradient
            else field_factory.uniform(strength)
        )
        builder = HamiltonianBuilder(
            protein=protein,
            interaction=interaction,
            distance_map=distance_map,
            contact_map=contact_map,
            external_field=field,
        )
        hamiltonian: SparsePauliOp = builder.sum_hamiltonians()
        compressed, kept = compress_with_layout(hamiltonian)
        exact: ExactSolution = solve_exactly(compressed)

        structure: DecodedStructure = decode_structure(
            bitstring=expand_bitstring(
                exact.bitstring, kept, int(hamiltonian.num_qubits)
            ),
            chain_length=len(sequence),
        )
        projections: NDArray[np.float64] = structure.coordinates @ field.direction
        extent: float = float(np.ptp(projections))

        records.append(
            {
                "mode": field.mode.value,
                "strength": strength,
                "energy": exact.energy,
                "degeneracy": exact.degeneracy,
                "turns": " ".join(str(turn.value) for turn in structure.turns),
                "extent_along_field": extent,
            }
        )
        logger.info(
            "Field %s strength %.3f -> E=%.4f, turns=%s",
            field.mode.value,
            strength,
            exact.energy,
            [turn.value for turn in structure.turns],
        )

    return records


def sweep_particle_environment(
    sequence: str,
    energy_scales: list[float],
    potential: RadialPotential | None = None,
    num_ligand_steps: int = DEFAULT_LIGAND_STEPS,
    hp_type: str = "H",
) -> list[SweepRecord]:
    """Trace how the chain arranges itself around a full-bead particle.

    Sweeping the interaction energy through zero covers both regimes in one
    pass. Note that *energy_scale* multiplies an already-negative contact energy,
    so a positive scale is an attraction and a negative one a repulsion; each
    record carries the resolved per-residue energy so the sign is never in doubt.
    Nothing in the Hamiltonian rewards coordination directly, so the coordination
    number reported here is emergent.

    Args:
        sequence (str): Main chain sequence.
        energy_scales (list[float]): Affinity multipliers applied to the contact
            energy. Negative values invert the interaction into a repulsion.
        potential (RadialPotential | None, optional): Shell profile. Defaults to
            the standard finite-range well.
        num_ligand_steps (int, optional): Lattice steps encoding the particle's
            position. Defaults to DEFAULT_LIGAND_STEPS.
        hp_type (str, optional): The particle's HP character. Defaults to "H".

    Returns:
        list[SweepRecord]: One record per affinity.

    """
    records: list[SweepRecord] = []
    side_chain: str = EMPTY_SIDECHAIN_PLACEHOLDER * len(sequence)
    protein, interaction, contact_map, distance_map = setup_folding_system(
        main_chain=sequence, side_chain=side_chain
    )
    ligand = Ligand(num_steps=num_ligand_steps)
    shell_profile: RadialPotential = potential or RadialPotential()

    for scale in energy_scales:
        ligand_interaction = LigandInteraction.hp_like(
            hp_type=hp_type, energy_scale=scale
        )
        strongest_energy: float = max(
            (ligand_interaction.get_energy(symbol) for symbol in sequence),
            key=abs,
        )
        builder = HamiltonianBuilder(
            protein=protein,
            interaction=interaction,
            distance_map=distance_map,
            contact_map=contact_map,
            ligand=ligand,
            ligand_interaction=ligand_interaction,
            coupling=LigandCoupling.RADIAL_POTENTIAL,
            radial_potential=shell_profile,
        )
        hamiltonian: SparsePauliOp = builder.sum_hamiltonians()
        compressed, kept = compress_with_layout(hamiltonian)
        exact: ExactSolution = solve_exactly(compressed)

        structure: DecodedStructure = decode_structure(
            bitstring=expand_bitstring(
                exact.bitstring, kept, int(hamiltonian.num_qubits)
            ),
            chain_length=len(sequence),
            ligand=ligand,
            ligand_walk_offset=builder.ligand_walk_offset,
            ligand_contact_offset=builder.contact_flag_offset,
        )
        environment: ParticleEnvironment = analyse_environment(
            hamiltonian=hamiltonian,
            chain_length=len(sequence),
            ligand=ligand,
            walk_qubit_offset=builder.ligand_walk_offset,
        )

        separations = structure.ligand_separations()

        records.append(
            {
                "energy_scale": scale,
                "interaction_energy": strongest_energy,
                "regime": (
                    "attractive"
                    if strongest_energy < 0
                    else ("repulsive" if strongest_energy > 0 else "neutral")
                ),
                "energy": exact.energy,
                "degeneracy": exact.degeneracy,
                "coordination_number": environment.coordination_number,
                "buried_fraction": environment.buried_fraction,
                "mean_coordination": environment.mean_coordination,
                "position_entropy": environment.position_entropy,
                "nearest_distance": float(separations.min()),
                "mean_distance": float(separations.mean()),
                "encapsulated": environment.is_encapsulated,
                "expelled": environment.is_expelled,
                "turns": " ".join(str(turn.value) for turn in structure.turns),
                "touching_residues": " ".join(
                    sequence[index] for index in structure.lattice_neighbours
                ),
            }
        )
        logger.info(
            "Affinity %.2f -> coordination %.0f, nearest %.3f, entropy %.3f",
            scale,
            environment.coordination_number,
            separations.min(),
            environment.position_entropy,
        )

    return records


def write_records_csv(records: list[SweepRecord], filepath: Path) -> Path:
    """Write sweep records to a CSV file.

    Args:
        records (list[SweepRecord]): Records sharing the same keys.
        filepath (Path): Destination file.

    Returns:
        Path: The path written to.

    Raises:
        ValueError: If *records* is empty.

    """
    if not records:
        msg: str = "No records to write"
        raise ValueError(msg)

    filepath.parent.mkdir(parents=True, exist_ok=True)
    with filepath.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    logger.info("Wrote %d records to %s", len(records), filepath)
    return filepath
