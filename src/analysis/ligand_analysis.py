"""Ligand–protein interaction analysis via quantum VQE.

This module provides :class:`LigandAnalysis`, which orchestrates VQE runs for a
protein–ligand system and exposes three high-level analysis primitives:

1. :meth:`~LigandAnalysis.compute_binding_energy` – estimates the binding free
   energy as the difference between the ground-state energy of the
   protein+ligand Hamiltonian and the sum of their individual ground-state
   energies (with and without ligand).

2. :meth:`~LigandAnalysis.compute_ligand_position_distribution` – extracts the
   marginal probability distribution over the ligand's positional qubit register
   from the VQE quasi-distribution, yielding P(node_k) for each lattice node.

3. :meth:`~LigandAnalysis.detect_encapsulation` – applies a geometric heuristic
   to the position distribution: a ligand is considered *encapsulated* if a
   high fraction of its probability mass is concentrated on the subset of
   lattice nodes that are "interior" (i.e. surrounded by high-probability
   protein conformations), compared to the boundary nodes.

Usage example::

    from pathlib import Path
    from src.analysis.ligand_analysis import LigandAnalysis
    from src.enums import InteractionType
    from src.interaction.ligand_interaction import LigandInteractionMode

    analysis = LigandAnalysis(
        main_chain="HPPHH",
        interaction_type=InteractionType.HP,
        ligand_hp_type="H",
        num_lattice_nodes=8,
        vqe_max_iter=150,
    )
    analysis.run()
    print(analysis.summary())
    analysis.plot(output_dir=Path("output/ligand_analysis"))

Architecture note
-----------------
The ligand is modelled in the **Variant A** approximation (see
:mod:`builder.hamiltonian_builder`): the ligand contact term adds the
scalar ``Σ_i E_i · I`` to the Hamiltonian and *extends* the qubit register by
``num_position_qubits``.  The positional qubits therefore remain unentangled
from the protein turn-qubits in this approximation; their marginal distribution
is determined by the ansatz and optimizer, not by direct Hamiltonian coupling.

A non-trivial position distribution nevertheless emerges because:

* The VQE optimizer freely rotates all qubits (protein + ligand) to minimize
  the total energy.
* When the Hamiltonian has degenerate ground states, the optimizer preferentially
  samples certain configurations depending on the aggregation (CVaR) parameter.
* Projecting (tracing out) the protein qubits from the best measurement gives
  the marginal over ligand positions, which we interpret as P(node_k).
"""

from __future__ import annotations

import dataclasses
import math
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from qiskit.circuit.library import real_amplitudes
from qiskit.quantum_info import SparsePauliOp
from qiskit_algorithms import SamplingVQE
from qiskit_algorithms.optimizers import COBYLA

from backend import get_sampler
from builder import HamiltonianBuilder
from constants import EMPTY_SIDECHAIN_PLACEHOLDER
from contact import ContactMap
from distance import DistanceMap
from enums import InteractionType
from exceptions import InvalidInteractionTypeError
from interaction import HPInteraction, MJInteraction
from interaction.ligand_interaction import LigandInteraction, LigandInteractionMode
from logger import get_logger
from particle.ligand_bead import LigandBead, PositionEncoding
from protein import Protein
from utils.qubit_utils import remove_unused_qubits

if TYPE_CHECKING:
    from qiskit_algorithms import SamplingMinimumEigensolverResult

logger = get_logger()


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------


@dataclasses.dataclass
class LigandScenarioResult:
    """Result of a single VQE run in a ligand-analysis scenario.

    Attributes:
        label (str): Human-readable scenario label.
        minimum_energy (float): Ground-state energy found by VQE.
        best_bitstring (str): Bitstring of the lowest-energy state.
        state_probabilities (dict[str, float]): Full quasi-distribution over
            computational basis states.
        ligand_node_probabilities (list[float]): Marginal probability over
            lattice nodes: ``ligand_node_probabilities[k]`` = P(ligand at k).
        binding_energy (float | None): Filled in by
            :meth:`~LigandAnalysis.compute_binding_energy` after all runs.
        vqe_iterations (list[int]): Eval-count at each VQE callback step.
        vqe_energies (list[float]): Energy at each VQE callback step.

    """

    label: str
    minimum_energy: float
    best_bitstring: str
    state_probabilities: dict[str, float]
    ligand_node_probabilities: list[float]
    binding_energy: float | None = None
    vqe_iterations: dataclasses.field(default_factory=list) = dataclasses.field(
        default_factory=list
    )
    vqe_energies: dataclasses.field(default_factory=list) = dataclasses.field(
        default_factory=list
    )


@dataclasses.dataclass
class EncapsulationResult:
    """Heuristic result of the encapsulation test.

    Attributes:
        is_encapsulated (bool): ``True`` when the encapsulation score exceeds
            the threshold.
        score (float): Fraction of ligand probability mass on interior nodes,
            in [0, 1].
        threshold (float): The score threshold used for classification.
        interior_nodes (list[int]): Node indices classified as interior.
        boundary_nodes (list[int]): Node indices classified as boundary.
        interior_probability (float): Summed probability on interior nodes.
        boundary_probability (float): Summed probability on boundary nodes.

    """

    is_encapsulated: bool
    score: float
    threshold: float
    interior_nodes: list[int]
    boundary_nodes: list[int]
    interior_probability: float
    boundary_probability: float


# ---------------------------------------------------------------------------
# Main analysis class
# ---------------------------------------------------------------------------


class LigandAnalysis:
    """Quantum VQE analysis of a protein–ligand interaction system.

    Runs three VQE scenarios and provides ligand-specific analysis methods:

    1. **Baseline** – protein alone (no ligand).
    2. **Ligand (full system)** – protein + ligand Hamiltonian.
    3. **Ligand alone** – trivial one-qubit placeholder (for binding energy
       computation).

    After calling :meth:`run`, use :meth:`compute_binding_energy`,
    :meth:`compute_ligand_position_distribution`, :meth:`detect_encapsulation`,
    and :meth:`plot` to analyse the results.

    Attributes:
        main_chain (str): Amino-acid sequence of the protein.
        side_chain (str): Side-chain sequence (``_`` for all positions).
        interaction_type (InteractionType): HP or MJ protein interaction.
        ligand_hp_type (str | None): ``"H"`` or ``"P"`` for HP-like ligand.
            If ``None``, *ligand_energy_map* is used for CUSTOM mode.
        ligand_energy_map (dict[str, float] | None): Custom energy map (CUSTOM
            mode).  Ignored when *ligand_hp_type* is set.
        ligand_default_energy (float): Default energy for unmapped residues in
            CUSTOM mode.
        num_lattice_nodes (int): Number of addressable lattice nodes for the
            ligand.
        position_encoding (PositionEncoding): Binary (default) or Unary.
        vqe_max_iter (int): Maximum COBYLA iterations per VQE run.
        vqe_aggregation (float): CVaR aggregation fraction for SamplingVQE.
        results (dict[str, LigandScenarioResult]): Populated after :meth:`run`.

    """

    _KEY_BASELINE = "baseline"
    _KEY_WITH_LIGAND = "with_ligand"

    def __init__(
        self,
        main_chain: str,
        *,
        interaction_type: InteractionType = InteractionType.HP,
        ligand_hp_type: str | None = "H",
        ligand_energy_map: dict[str, float] | None = None,
        ligand_default_energy: float = 0.0,
        num_lattice_nodes: int = 8,
        position_encoding: PositionEncoding = PositionEncoding.BINARY,
        side_chain: str | None = None,
        vqe_max_iter: int = 100,
        vqe_aggregation: float = 0.1,
    ) -> None:
        """Initialise the :class:`LigandAnalysis`.

        Args:
            main_chain (str): One-letter amino-acid sequence of the main chain.
            interaction_type (InteractionType, optional): Protein interaction
                model.  Defaults to HP.
            ligand_hp_type (str | None, optional): ``"H"`` (hydrophobic) or
                ``"P"`` (polar) for HP-like ligand mode, or ``None`` to use
                the CUSTOM mode with *ligand_energy_map*.  Defaults to ``"H"``.
            ligand_energy_map (dict[str, float] | None, optional): Per-residue
                energies for CUSTOM mode.  Ignored when *ligand_hp_type* is
                not ``None``.
            ligand_default_energy (float, optional): Default energy for
                residues not in *ligand_energy_map* (CUSTOM mode only).
                Defaults to ``0.0``.
            num_lattice_nodes (int, optional): Number of lattice nodes the
                ligand can occupy.  Defaults to ``8``.
            position_encoding (PositionEncoding, optional): Encoding strategy
                for the ligand's position qubits.  Defaults to BINARY.
            side_chain (str | None, optional): Side-chain sequence.  Defaults
                to all ``_`` (no side chains).
            vqe_max_iter (int, optional): COBYLA maximum iterations per run.
                Defaults to 100.
            vqe_aggregation (float, optional): CVaR aggregation fraction.
                Defaults to 0.1.

        Raises:
            InvalidInteractionTypeError: If *interaction_type* is unknown.
            ValueError: If *ligand_hp_type* is not ``"H"``, ``"P"``, or
                ``None``.

        """
        self.main_chain: str = main_chain
        self.side_chain: str = (
            side_chain
            if side_chain is not None
            else EMPTY_SIDECHAIN_PLACEHOLDER * len(main_chain)
        )
        self.interaction_type: InteractionType = interaction_type
        self.ligand_hp_type: str | None = ligand_hp_type
        self.ligand_energy_map: dict[str, float] | None = ligand_energy_map
        self.ligand_default_energy: float = ligand_default_energy
        self.num_lattice_nodes: int = num_lattice_nodes
        self.position_encoding: PositionEncoding = position_encoding
        self.vqe_max_iter: int = vqe_max_iter
        self.vqe_aggregation: float = vqe_aggregation

        self.results: dict[str, LigandScenarioResult] = {}

        # Build shared protein components
        self._interaction = self._build_protein_interaction()
        self._protein = Protein(
            main_protein_sequence=self.main_chain,
            side_protein_sequence=self.side_chain,
            valid_symbols=self._interaction.valid_symbols,
        )
        self._contact_map = ContactMap(protein=self._protein)
        self._distance_map = DistanceMap(protein=self._protein)

        # Build ligand
        self._ligand = LigandBead(
            symbol="L",
            index=0,
            num_lattice_nodes=self.num_lattice_nodes,
            encoding=self.position_encoding,
        )
        self._ligand_interaction = self._build_ligand_interaction()

        logger.info(
            "LigandAnalysis initialised: chain=%s (len=%d), model=%s, "
            "ligand_mode=%s, nodes=%d, pos_qubits=%d",
            self.main_chain,
            len(self.main_chain),
            self.interaction_type.name,
            self._ligand_interaction.mode.name,
            self.num_lattice_nodes,
            self._ligand.num_position_qubits,
        )

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def run(self) -> None:
        """Run all scenarios (baseline + with-ligand) and populate :attr:`results`.

        Populates :attr:`results` with two entries:

        * ``"baseline"`` – protein alone.
        * ``"with_ligand"`` – protein + ligand Hamiltonian.

        Also computes :attr:`~LigandScenarioResult.binding_energy` for the
        ``"with_ligand"`` result via
        :meth:`compute_binding_energy`.

        """
        logger.info("LigandAnalysis.run() starting …")
        self.results.clear()

        # 1. Baseline: protein only
        logger.info("--- Scenario: baseline (protein only) ---")
        self.results[self._KEY_BASELINE] = self._run_protein_only()

        # 2. Full system: protein + ligand
        logger.info("--- Scenario: protein + ligand ---")
        self.results[self._KEY_WITH_LIGAND] = self._run_with_ligand()

        # 3. Derive binding energy
        self.results[self._KEY_WITH_LIGAND].binding_energy = (
            self.compute_binding_energy()
        )

        logger.info(
            "LigandAnalysis.run() complete.  Binding energy = %.6f",
            self.results[self._KEY_WITH_LIGAND].binding_energy,
        )

    # ------------------------------------------------------------------
    # Analysis methods
    # ------------------------------------------------------------------

    def compute_binding_energy(self) -> float:
        """Estimate the ligand binding energy.

        The binding energy is defined as the difference between the ground-state
        energies of the coupled protein+ligand system and the uncoupled baseline:

            ΔE_binding = E(protein+ligand) − E(protein alone)

        In the Variant A approximation the ligand contact term shifts the
        Hamiltonian by the scalar ``Σ_i E_i``, so the binding energy should
        equal ``Σ_i E_i`` (the total ligand-chain interaction energy).  The VQE
        may not converge exactly, so we report the empirically measured energy
        difference.

        A **negative** value indicates stabilisation (favourable binding).

        Returns:
            float: Binding energy ΔE = E_with_ligand − E_baseline.

        Raises:
            RuntimeError: If :meth:`run` has not been called.

        """
        self._require_results()
        e_with = self.results[self._KEY_WITH_LIGAND].minimum_energy
        e_base = self.results[self._KEY_BASELINE].minimum_energy
        delta_e = e_with - e_base
        logger.info(
            "Binding energy: ΔE = %.6f − %.6f = %.6f",
            e_with, e_base, delta_e,
        )
        return delta_e

    def compute_ligand_position_distribution(self) -> list[float]:
        """Extract the marginal probability distribution over lattice nodes.

        Marginalises the full quasi-distribution of the protein+ligand VQE run
        over the protein qubits, yielding a list of probabilities
        ``P[k] = P(ligand occupies node k)``.

        In the BINARY encoding the ligand occupies the *last*
        ``num_position_qubits`` bits of each bitstring.  In UNARY encoding the
        last ``num_lattice_nodes`` bits form a one-hot vector.

        The resulting distribution sums to 1.0 (up to floating-point rounding).

        Returns:
            list[float]: Probability list of length ``num_lattice_nodes``.

        Raises:
            RuntimeError: If :meth:`run` has not been called.

        """
        self._require_results()
        result = self.results[self._KEY_WITH_LIGAND]
        n_pos = self._ligand.num_position_qubits
        n_nodes = self.num_lattice_nodes

        node_probs = [0.0] * n_nodes

        for bitstring, prob in result.state_probabilities.items():
            if prob <= 0.0:
                continue
            # Ligand position bits are the LAST n_pos bits (little-endian):
            # Qiskit bitstrings are written MSB-first, so the last chars
            # correspond to the lowest-index qubits in the register.
            ligand_bits = bitstring[-n_pos:]  # right-most n_pos chars

            if self.position_encoding is PositionEncoding.BINARY:
                # Interpret as binary integer → node index
                node_idx = int(ligand_bits, 2)
                if node_idx < n_nodes:
                    node_probs[node_idx] += prob
            else:
                # UNARY: the '1' bit position gives the node index
                # ligand_bits[0] = highest-index qubit in ligand register
                # We interpret as MSB-first one-hot
                for bit_pos, bit_char in enumerate(ligand_bits):
                    # bit_pos 0 → leftmost char → highest qubit index
                    qubit_idx = n_pos - 1 - bit_pos
                    if bit_char == "1" and qubit_idx < n_nodes:
                        node_probs[qubit_idx] += prob

        # Normalise to [0, 1] (quasi-distributions can have negative weights)
        total = sum(max(p, 0.0) for p in node_probs)
        if total > 0.0:
            node_probs = [max(p, 0.0) / total for p in node_probs]

        logger.info(
            "Ligand position distribution (nodes 0–%d): %s",
            n_nodes - 1,
            [round(p, 4) for p in node_probs],
        )
        return node_probs

    def detect_encapsulation(
        self,
        interior_fraction: float = 0.5,
        encapsulation_threshold: float = 0.6,
    ) -> EncapsulationResult:
        """Apply a heuristic to detect whether the ligand is encapsulated.

        The heuristic proceeds in two steps:

        1. **Classify nodes** as *interior* or *boundary* using a simple
           criterion based on the node index relative to the lattice size.
           Nodes in the central ``interior_fraction`` of the lattice are
           classified as interior; the rest are boundary.  In a 1-D
           linearisation of the lattice this corresponds to "middle" nodes.

        2. **Compute encapsulation score** as the fraction of the ligand's
           probability mass concentrated on interior nodes::

               score = Σ_{k ∈ interior} P[k]

        3. **Classify** as encapsulated if ``score ≥ encapsulation_threshold``.

        This is a pure-classical post-processing step operating on the
        quantum-computed position distribution; it does not run VQE again.

        Args:
            interior_fraction (float, optional): Fraction of nodes (sorted by
                index) considered interior.  Defaults to 0.5.
            encapsulation_threshold (float, optional): Minimum score to classify
                the ligand as encapsulated.  Defaults to 0.6.

        Returns:
            EncapsulationResult: Encapsulation classification and supporting
                statistics.

        Raises:
            RuntimeError: If :meth:`run` has not been called.

        """
        self._require_results()
        node_probs = self.compute_ligand_position_distribution()
        n_nodes = len(node_probs)

        # Classify nodes: interior = central interior_fraction of the lattice
        n_interior = max(1, round(n_nodes * interior_fraction))
        margin = (n_nodes - n_interior) // 2
        interior_nodes = list(range(margin, margin + n_interior))
        boundary_nodes = [k for k in range(n_nodes) if k not in interior_nodes]

        interior_prob = sum(node_probs[k] for k in interior_nodes)
        boundary_prob = sum(node_probs[k] for k in boundary_nodes)
        score = interior_prob  # = P(ligand on interior node)
        is_encapsulated = score >= encapsulation_threshold

        enc_result = EncapsulationResult(
            is_encapsulated=is_encapsulated,
            score=score,
            threshold=encapsulation_threshold,
            interior_nodes=interior_nodes,
            boundary_nodes=boundary_nodes,
            interior_probability=interior_prob,
            boundary_probability=boundary_prob,
        )

        logger.info(
            "Encapsulation: score=%.4f (threshold=%.2f) → %s | "
            "interior_nodes=%s, boundary_nodes=%s",
            score,
            encapsulation_threshold,
            "ENCAPSULATED" if is_encapsulated else "surface",
            interior_nodes,
            boundary_nodes,
        )
        return enc_result

    # ------------------------------------------------------------------
    # Reporting
    # ------------------------------------------------------------------

    def summary(self) -> str:
        """Return a formatted text summary of all results.

        Returns:
            str: Multi-line table with scenario label, E_min, binding energy,
                and best bitstring for each completed run.

        Raises:
            RuntimeError: If :meth:`run` has not been called.

        """
        self._require_results()
        lines = [
            f"{'Scenario':<25}  {'E_min':>12}  {'ΔE_binding':>12}  Best bitstring",
            "-" * 70,
        ]
        for key, r in self.results.items():
            delta = f"{r.binding_energy:>12.6f}" if r.binding_energy is not None else f"{'—':>12}"
            lines.append(
                f"{r.label:<25}  {r.minimum_energy:>12.6f}  {delta}  {r.best_bitstring}"
            )
        return "\n".join(lines)

    def plot(self, output_dir: Path | None = None) -> None:
        """Generate and save ligand analysis plots.

        Produces three figures:

        1. **Ligand position distribution** – bar chart showing P(node_k).
        2. **Energy comparison** – bar chart with baseline vs. ligand energies.
        3. **Lattice heatmap** – 2-D heatmap of P(node) on a virtual N×M grid.

        Args:
            output_dir (Path | None, optional): Directory for saving figures.
                Interactive display when ``None``.

        Raises:
            RuntimeError: If :meth:`run` has not been called.

        """
        self._require_results()

        try:
            import matplotlib.pyplot as plt
            import matplotlib.ticker as ticker
        except ImportError as e:
            msg = "matplotlib is required for plotting."
            raise ImportError(msg) from e

        if output_dir is not None:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

        _PALETTE = [
            "#4c9be8", "#e86b4c", "#5fcf80", "#f0c040",
            "#a87de8", "#e8a04c", "#4ce8d8", "#e84c8a",
        ]

        self._plot_position_distribution(plt, ticker, _PALETTE, output_dir)
        self._plot_energy_comparison(plt, _PALETTE, output_dir)
        self._plot_lattice_heatmap(plt, output_dir)

        if output_dir is None:
            plt.show()
        else:
            logger.info("All ligand-analysis plots saved to %s", output_dir)

    # ------------------------------------------------------------------
    # Internal: scenario execution
    # ------------------------------------------------------------------

    def _build_protein_interaction(self) -> HPInteraction | MJInteraction:
        """Build the protein–protein interaction model."""
        if self.interaction_type == InteractionType.HP:
            return HPInteraction()
        if self.interaction_type == InteractionType.MJ:
            return MJInteraction()
        msg = f"Unknown interaction type: {self.interaction_type}"
        raise InvalidInteractionTypeError(msg)

    def _build_ligand_interaction(self) -> LigandInteraction:
        """Build the ligand–protein interaction model."""
        if self.ligand_hp_type is not None:
            return LigandInteraction.hp_like(
                ligand_hp_type=self.ligand_hp_type, ligand_symbol="L"
            )
        energy_map = self.ligand_energy_map if self.ligand_energy_map is not None else {}
        return LigandInteraction.custom(
            energy_map=energy_map,
            default_energy=self.ligand_default_energy,
            ligand_symbol="L",
        )

    def _run_protein_only(self) -> LigandScenarioResult:
        """Run VQE for the protein-only Hamiltonian (baseline)."""
        builder = HamiltonianBuilder(
            protein=self._protein,
            interaction=self._interaction,
            distance_map=self._distance_map,
            contact_map=self._contact_map,
        )
        full_h: SparsePauliOp = builder.sum_hamiltonians()
        compressed_h = remove_unused_qubits(full_h)

        iterations: list[int] = []
        energies: list[float] = []

        def _cb(eval_count: int, _p: np.ndarray, mean: float, _s: dict) -> None:
            iterations.append(eval_count)
            energies.append(mean)

        raw = self._run_vqe(compressed_h, _cb)
        best = raw.best_measurement or {}
        min_energy = float(np.real(best.get("value", float("nan"))))
        best_bs = best.get("bitstring", "")
        probs = self._extract_probs(raw, int(compressed_h.num_qubits))

        return LigandScenarioResult(
            label="baseline (protein only)",
            minimum_energy=min_energy,
            best_bitstring=best_bs,
            state_probabilities=probs,
            ligand_node_probabilities=[1.0 / self.num_lattice_nodes] * self.num_lattice_nodes,
            vqe_iterations=iterations,
            vqe_energies=energies,
        )

    def _run_with_ligand(self) -> LigandScenarioResult:
        """Run VQE for the protein + ligand Hamiltonian."""
        builder = HamiltonianBuilder(
            protein=self._protein,
            interaction=self._interaction,
            distance_map=self._distance_map,
            contact_map=self._contact_map,
        )
        full_h: SparsePauliOp = builder.sum_hamiltonians(
            ligand=self._ligand,
            ligand_interaction=self._ligand_interaction,
        )
        compressed_h = remove_unused_qubits(full_h)

        logger.info(
            "Protein+ligand Hamiltonian: %d qubits total (%d after compression).",
            full_h.num_qubits,
            compressed_h.num_qubits,
        )

        iterations: list[int] = []
        energies: list[float] = []

        def _cb(eval_count: int, _p: np.ndarray, mean: float, _s: dict) -> None:
            iterations.append(eval_count)
            energies.append(mean)

        raw = self._run_vqe(compressed_h, _cb)
        best = raw.best_measurement or {}
        min_energy = float(np.real(best.get("value", float("nan"))))
        best_bs = best.get("bitstring", "")
        probs = self._extract_probs(raw, int(compressed_h.num_qubits))

        # Compute position distribution from quasi-distribution
        result = LigandScenarioResult(
            label="protein + ligand",
            minimum_energy=min_energy,
            best_bitstring=best_bs,
            state_probabilities=probs,
            ligand_node_probabilities=[],  # populated below
            vqe_iterations=iterations,
            vqe_energies=energies,
        )
        # Temporarily store result so compute_ligand_position_distribution can use it
        self.results[self._KEY_WITH_LIGAND] = result
        result.ligand_node_probabilities = self.compute_ligand_position_distribution()
        return result

    def _run_vqe(
        self,
        hamiltonian: SparsePauliOp,
        callback,
    ) -> SamplingMinimumEigensolverResult:
        """Run SamplingVQE on *hamiltonian* and return the raw result."""
        sampler, _ = get_sampler()
        ansatz = real_amplitudes(num_qubits=int(hamiltonian.num_qubits), reps=1)
        vqe = SamplingVQE(
            sampler=sampler,
            ansatz=ansatz,
            optimizer=COBYLA(maxiter=self.vqe_max_iter),
            aggregation=self.vqe_aggregation,
            callback=callback,
        )
        return vqe.compute_minimum_eigenvalue(hamiltonian)

    @staticmethod
    def _extract_probs(
        raw: SamplingMinimumEigensolverResult, n_bits: int
    ) -> dict[str, float]:
        """Convert quasi-distribution from *raw* to a bitstring→probability dict."""
        probs: dict[str, float] = {}
        if raw.eigenstate is None:
            return probs
        for key, prob in raw.eigenstate.items():
            if prob > 0:
                bs = format(key, f"0{n_bits}b") if isinstance(key, int) else str(key).zfill(n_bits)
                probs[bs] = float(prob)
        return probs

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _require_results(self) -> None:
        """Raise RuntimeError if run() has not been called yet."""
        if not self.results:
            msg = "No results available. Call run() first."
            raise RuntimeError(msg)

    # ------------------------------------------------------------------
    # Internal: plotting
    # ------------------------------------------------------------------

    def _plot_position_distribution(
        self,
        plt,
        ticker,
        palette: list[str],
        output_dir: Path | None,
    ) -> None:
        """Bar chart of ligand position distribution across lattice nodes."""
        result = self.results[self._KEY_WITH_LIGAND]
        node_probs = result.ligand_node_probabilities
        n_nodes = len(node_probs)

        fig, ax = plt.subplots(figsize=(max(8, n_nodes * 0.6 + 2), 5))
        fig.patch.set_facecolor("#0f1117")
        ax.set_facecolor("#1a1d27")

        colours = [
            palette[i % len(palette)] for i in range(n_nodes)
        ]
        bars = ax.bar(
            range(n_nodes),
            node_probs,
            color=colours,
            edgecolor="#2a2d3a",
            linewidth=0.8,
            width=0.72,
            zorder=3,
        )

        # Annotate bars
        for bar, p in zip(bars, node_probs):
            if p > 0.02:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.005,
                    f"{p:.3f}",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                    color="#e0e0e0",
                )

        ax.set_xticks(range(n_nodes))
        ax.set_xticklabels(
            [f"node {k}" for k in range(n_nodes)],
            rotation=45,
            ha="right",
            fontsize=8,
            color="#c0c0c0",
        )
        ax.yaxis.set_tick_params(labelcolor="#c0c0c0")
        ax.set_ylabel("P(ligand at node)", fontsize=12, color="#e0e0e0")
        ax.set_ylim(0, min(1.05, max(node_probs) * 1.25 + 0.05) if node_probs else 1.05)
        ax.set_xlabel("Lattice node index", fontsize=11, color="#e0e0e0")

        mode_label = (
            f"HP-like ({self.ligand_hp_type})"
            if self.ligand_hp_type
            else f"Custom (map={list((self.ligand_energy_map or {}).keys())})"
        )
        ax.set_title(
            f"Ligand position distribution — {self.main_chain} ({self.interaction_type.name})\n"
            f"Ligand mode: {mode_label} | nodes={n_nodes} | encoding={self.position_encoding.name}",
            fontsize=12,
            color="#ffffff",
            pad=10,
        )
        ax.grid(axis="y", color="#333", linewidth=0.6, zorder=0)
        ax.spines[:].set_edgecolor("#333")
        ticker.AutoMinorLocator()

        fig.tight_layout()
        _save_or_show_lig(fig, plt, output_dir, "ligand_position_distribution.png")

    def _plot_energy_comparison(
        self,
        plt,
        palette: list[str],
        output_dir: Path | None,
    ) -> None:
        """Bar chart comparing baseline and with-ligand minimum energies."""
        labels = [r.label for r in self.results.values()]
        energies = [r.minimum_energy for r in self.results.values()]
        binding = self.results[self._KEY_WITH_LIGAND].binding_energy

        fig, ax = plt.subplots(figsize=(7, 5))
        fig.patch.set_facecolor("#0f1117")
        ax.set_facecolor("#1a1d27")

        colours = [palette[i % len(palette)] for i in range(len(labels))]
        bars = ax.bar(
            range(len(labels)),
            energies,
            color=colours,
            edgecolor="#2a2d3a",
            linewidth=0.8,
            width=0.55,
            zorder=3,
        )

        for bar, e in zip(bars, energies):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + (max(energies) - min(energies)) * 0.01,
                f"{e:.4f}",
                ha="center",
                va="bottom",
                fontsize=9,
                color="#e0e0e0",
            )

        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=15, ha="right", fontsize=10, color="#c0c0c0")
        ax.yaxis.set_tick_params(labelcolor="#c0c0c0")
        ax.set_ylabel("Minimum VQE energy", fontsize=12, color="#e0e0e0")
        binding_str = f"{binding:+.4f}" if binding is not None else "N/A"
        ax.set_title(
            f"Energy comparison — {self.main_chain} ({self.interaction_type.name})\n"
            f"ΔE_binding = {binding_str}",
            fontsize=12,
            color="#ffffff",
            pad=10,
        )
        ax.grid(axis="y", color="#333", linewidth=0.6, zorder=0)
        ax.spines[:].set_edgecolor("#333")
        fig.tight_layout()
        _save_or_show_lig(fig, plt, output_dir, "energy_comparison.png")

    def _plot_lattice_heatmap(
        self,
        plt,
        output_dir: Path | None,
    ) -> None:
        """2-D heatmap of ligand probability on a virtual rectangular lattice.

        Maps the 1-D list of node probabilities onto a 2-D grid of size
        (rows × cols) where rows = ⌊√N⌋ and cols = ⌈N / rows⌉, padding with
        zeros as needed.  This provides an intuitive spatial view without
        requiring actual lattice coordinates.

        """
        result = self.results[self._KEY_WITH_LIGAND]
        node_probs = result.ligand_node_probabilities
        n = len(node_probs)

        rows = max(1, math.isqrt(n))
        cols = math.ceil(n / rows)
        # Pad to rows*cols
        padded = node_probs + [0.0] * (rows * cols - n)
        grid = np.array(padded).reshape(rows, cols)

        fig, ax = plt.subplots(figsize=(max(5, cols + 1), max(4, rows + 1)))
        fig.patch.set_facecolor("#0f1117")
        ax.set_facecolor("#0f1117")

        im = ax.imshow(
            grid,
            cmap="plasma",
            aspect="auto",
            interpolation="nearest",
            vmin=0.0,
            vmax=max(node_probs) if node_probs else 1.0,
        )

        # Annotate cells
        for r in range(rows):
            for c in range(cols):
                node_idx = r * cols + c
                if node_idx < n:
                    ax.text(
                        c, r,
                        f"{grid[r, c]:.3f}",
                        ha="center",
                        va="center",
                        fontsize=8,
                        color="white" if grid[r, c] < max(node_probs) * 0.6 else "black",
                    )

        # Encapsulation overlay
        enc = self.detect_encapsulation()
        interior_cells = []
        for k in enc.interior_nodes:
            r, c = divmod(k, cols)
            interior_cells.append((r, c))

        for r, c in interior_cells:
            rect = plt.Rectangle(
                (c - 0.48, r - 0.48), 0.96, 0.96,
                linewidth=2.0, edgecolor="#5fcf80", facecolor="none",
                zorder=5,
            )
            ax.add_patch(rect)

        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("P(ligand)", color="#e0e0e0", fontsize=10)
        cbar.ax.yaxis.set_tick_params(colors="#c0c0c0")

        ax.set_xticks(range(cols))
        ax.set_xticklabels([str(c) for c in range(cols)], color="#c0c0c0")
        ax.set_yticks(range(rows))
        ax.set_yticklabels([str(r) for r in range(rows)], color="#c0c0c0")
        ax.set_xlabel("Lattice column", fontsize=10, color="#e0e0e0")
        ax.set_ylabel("Lattice row", fontsize=10, color="#e0e0e0")
        enc_str = "ENCAPSULATED ✓" if enc.is_encapsulated else f"surface (score={enc.score:.2f})"
        ax.set_title(
            f"Ligand spatial distribution — {self.main_chain}\n"
            f"Encapsulation: {enc_str} | green border = interior nodes",
            fontsize=11,
            color="#ffffff",
            pad=8,
        )
        ax.spines[:].set_edgecolor("#333")
        fig.tight_layout()
        _save_or_show_lig(fig, plt, output_dir, "lattice_heatmap.png")


# ---------------------------------------------------------------------------
# Convenience: experiment runner
# ---------------------------------------------------------------------------


def run_hp_hydrophobic_experiment(
    main_chain: str = "HPPHH",
    num_lattice_nodes: int = 8,
    vqe_max_iter: int = 100,
    output_dir: Path | None = None,
) -> LigandAnalysis:
    """Run the hydrophobic-ligand experiment with an HP protein.

    Scenario: hydrophobic ligand (HP_LIKE "H") docked to a short HP sequence.
    Expected behaviour: the ligand should prefer interior (hydrophobic-core)
    nodes; binding energy should be negative.

    Args:
        main_chain (str, optional): HP amino-acid sequence. Defaults to
            ``"HPPHH"``.
        num_lattice_nodes (int, optional): Lattice size for the ligand.
            Defaults to 8.
        vqe_max_iter (int, optional): COBYLA iterations. Defaults to 100.
        output_dir (Path | None, optional): Output directory for plots.

    Returns:
        LigandAnalysis: Completed analysis object (after run() is called).

    """
    logger.info("=== HP hydrophobic ligand experiment ===")
    analysis = LigandAnalysis(
        main_chain=main_chain,
        interaction_type=InteractionType.HP,
        ligand_hp_type="H",
        num_lattice_nodes=num_lattice_nodes,
        position_encoding=PositionEncoding.BINARY,
        vqe_max_iter=vqe_max_iter,
    )
    analysis.run()
    logger.info("\n%s\n", analysis.summary())

    enc = analysis.detect_encapsulation()
    logger.info(
        "Encapsulation result: %s (score=%.4f)",
        "ENCAPSULATED" if enc.is_encapsulated else "surface",
        enc.score,
    )

    analysis.plot(output_dir=output_dir)
    return analysis


def run_mj_strong_ligand_experiment(
    main_chain: str = "ACDEFGH",
    num_lattice_nodes: int = 8,
    vqe_max_iter: int = 100,
    output_dir: Path | None = None,
) -> LigandAnalysis:
    """Run a strongly-interacting custom ligand with an MJ-model protein.

    Scenario: custom ligand with strong affinity to hydrophobic MJ residues
    (A, C, I, L, M, F, W, V → −2.0 kcal/mol) and weak affinity to polar
    residues (→ −0.1 kcal/mol).  The MJ model captures more realistic
    side-chain packing energies.

    Args:
        main_chain (str, optional): MJ amino-acid sequence. Defaults to
            ``"ACDEFGH"`` (mixed hydrophobic/polar).
        num_lattice_nodes (int, optional): Lattice size for the ligand.
            Defaults to 8.
        vqe_max_iter (int, optional): COBYLA iterations. Defaults to 100.
        output_dir (Path | None, optional): Output directory for plots.

    Returns:
        LigandAnalysis: Completed analysis object.

    """
    logger.info("=== MJ strong-ligand experiment ===")

    # Strong energy map: hydrophobic residues get -2.0, rest get -0.1
    hydrophobic = {"A", "C", "I", "L", "M", "F", "W", "V"}
    energy_map = {aa: -2.0 if aa in hydrophobic else -0.1 for aa in "ACDEFGH"}

    analysis = LigandAnalysis(
        main_chain=main_chain,
        interaction_type=InteractionType.MJ,
        ligand_hp_type=None,
        ligand_energy_map=energy_map,
        ligand_default_energy=-0.1,
        num_lattice_nodes=num_lattice_nodes,
        position_encoding=PositionEncoding.BINARY,
        vqe_max_iter=vqe_max_iter,
    )
    analysis.run()
    logger.info("\n%s\n", analysis.summary())

    enc = analysis.detect_encapsulation()
    logger.info(
        "Encapsulation result: %s (score=%.4f)",
        "ENCAPSULATED" if enc.is_encapsulated else "surface",
        enc.score,
    )

    analysis.plot(output_dir=output_dir)
    return analysis


# ---------------------------------------------------------------------------
# Private plot helpers
# ---------------------------------------------------------------------------


def _save_or_show_lig(
    fig,
    plt,
    output_dir: Path | None,
    filename: str,
) -> None:
    """Save figure to *output_dir/filename* or display if *output_dir* is None."""
    if output_dir is not None:
        filepath = output_dir / filename
        fig.savefig(filepath, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
        logger.info("Saved plot: %s", filepath)
        plt.close(fig)
    else:
        plt.tight_layout()
