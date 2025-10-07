# (C) Copyright 2025 QFold-Thesis.
#
# This code is licensed under the Apache License, Version 2.0.
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""Visualization utilities for quantum protein folding results."""

from __future__ import annotations

from pathlib import Path
from typing import Optional
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from protein.protein import Protein
from constants import QUBITS_PER_TURN, EMPTY_SIDECHAIN_PLACEHOLDER
from logger import get_logger

logger = get_logger()


class ProteinFoldingVisualizer:
    """
    A comprehensive visualization tool for quantum protein folding results.
    
    This class handles:
    - Decoding VQE bitstring results into protein coordinates
    - Generating .xyz files for molecular visualization
    - Creating 3D plots with smooth tetrahedral lattice representation
    """

    # Tetrahedral lattice coordinates (4 directions from center)
    # Normalized coordinates of the 4 edges of a tetrahedron centered at 0
    TETRAHEDRAL_COORDINATES = (1.0 / np.sqrt(3)) * np.array([
        [-1, 1, 1],   # Direction 0
        [1, 1, -1],   # Direction 1
        [-1, -1, -1], # Direction 2
        [1, -1, 1]    # Direction 3
    ])

    def __init__(self, protein: Protein):
        """
        Initialize the visualizer with a protein structure.
        
        Args:
            protein: The protein object containing main and side chain information
        """
        self.protein = protein
        self.main_chain_sequence = protein.main_chain.sequence
        self.side_chain_sequence = protein.side_chain.sequence
        
        # Cache for computed coordinates
        self._main_positions = None
        self._side_positions = None
        
    def decode_vqe_result(self, raw_result: str, unused_qubits: Optional[List[int]] = None) -> Tuple[List[int], List[Optional[int]]]:
        """
        Decode VQE bitstring result into main chain and side chain turns.
        
        Args:
            raw_result: The bitstring result from VQE optimization
            unused_qubits: List of qubit indices that were unused during optimization
            
        Returns:
            Tuple of (main_chain_turns, side_chain_turns)
        """
        if unused_qubits is None:
            unused_qubits = []
            
        logger.info(f"Decoding VQE result: {raw_result}")
        
        # Convert bitstring to list of turns
        main_chain_turns = []
        side_chain_turns = []
        
        # Process bitstring in pairs (2 qubits per turn for DENSE encoding)
        bit_pairs = []
        for i in range(0, len(raw_result), QUBITS_PER_TURN):
            if i + 1 < len(raw_result):
                # Convert binary pair to decimal (0-3)
                turn_value = int(raw_result[i:i+2], 2)
                bit_pairs.append(turn_value)
        
        # Calculate number of main chain turns needed
        main_chain_length = len(self.main_chain_sequence)
        num_main_turns = main_chain_length - 1  # N beads need N-1 turns
        
        # Extract main chain turns
        main_chain_turns = bit_pairs[:num_main_turns] if bit_pairs else []
        
        # Extract side chain turns
        remaining_turns = bit_pairs[num_main_turns:] if len(bit_pairs) > num_main_turns else []
        
        # Map side chain turns to positions
        side_chain_turns = []
        turn_index = 0
        
        for i, side_char in enumerate(self.side_chain_sequence):
            if side_char != EMPTY_SIDECHAIN_PLACEHOLDER:
                if turn_index < len(remaining_turns):
                    side_chain_turns.append(remaining_turns[turn_index])
                    turn_index += 1
                else:
                    side_chain_turns.append(None)
            else:
                side_chain_turns.append(None)
                
        logger.info(f"Decoded main chain turns: {main_chain_turns}")
        logger.info(f"Decoded side chain turns: {side_chain_turns}")
        
        return main_chain_turns, side_chain_turns
    
    def generate_coordinates(self, main_chain_turns: List[int], side_chain_turns: List[Optional[int]]) -> Tuple[np.ndarray, List[Optional[np.ndarray]]]:
        """
        Generate 3D coordinates for main chain and side chains based on turns.
        
        Args:
            main_chain_turns: List of turn directions for main chain
            side_chain_turns: List of turn directions for side chains (None for no side chain)
            
        Returns:
            Tuple of (main_positions, side_positions)
        """
        logger.info("Generating 3D coordinates...")
        
        # Generate main chain positions
        main_positions = self._generate_main_chain_positions(main_chain_turns)
        
        # Generate side chain positions
        side_positions = self._generate_side_chain_positions(main_positions, side_chain_turns)
        
        self._main_positions = main_positions
        self._side_positions = side_positions
        
        return main_positions, side_positions
    
    def _generate_main_chain_positions(self, turns: List[int]) -> np.ndarray:
        """Generate main chain positions using tetrahedral lattice."""
        num_beads = len(self.main_chain_sequence)
        positions = np.zeros((num_beads, 3), dtype=float)
        
        # First bead at origin
        positions[0] = [0.0, 0.0, 0.0]
        
        # Generate subsequent positions based on turns
        for i, turn in enumerate(turns):
            if i + 1 < num_beads:
                # Alternate direction multiplier for tetrahedral lattice
                direction_multiplier = (-1) ** i
                direction_vector = direction_multiplier * self.TETRAHEDRAL_COORDINATES[turn]
                positions[i + 1] = positions[i] + direction_vector
                
        return positions
    
    def _generate_side_chain_positions(self, main_positions: np.ndarray, side_turns: List[Optional[int]]) -> List[Optional[np.ndarray]]:
        """Generate side chain positions."""
        side_positions = []
        
        for i, (main_pos, side_turn) in enumerate(zip(main_positions, side_turns)):
            if side_turn is None:
                side_positions.append(None)
            else:
                # Alternate direction for side chains
                direction_multiplier = (-1) ** (i + 1)
                side_vector = direction_multiplier * self.TETRAHEDRAL_COORDINATES[side_turn]
                side_position = main_pos + side_vector
                side_positions.append(side_position)
                
        return side_positions
    
    def save_xyz_file(self, 
                     filename: str, 
                     main_positions: Optional[np.ndarray] = None,
                     side_positions: Optional[List[Optional[np.ndarray]]] = None,
                     path: str = "", 
                     comment: str = "",
                     replace: bool = False) -> None:
        """
        Save protein coordinates as .xyz file.
        
        Args:
            filename: Name of the output file (without .xyz extension)
            main_positions: Main chain coordinates (uses cached if None)
            side_positions: Side chain coordinates (uses cached if None)
            path: Directory path to save file
            comment: Comment line for .xyz file
            replace: Whether to overwrite existing files
            
        Raises:
            FileExistsError: If file exists and replace is False
            ValueError: If no coordinates are available
        """
        if main_positions is None:
            main_positions = self._main_positions
        if side_positions is None:
            side_positions = self._side_positions
            
        if main_positions is None:
            raise ValueError("No coordinates available. Generate coordinates first.")
            
        file_path = os.path.join(path, f"{filename}.xyz")
        
        if not replace and os.path.exists(file_path):
            raise FileExistsError(f"File {file_path} already exists.")
            
        # Prepare data for .xyz format
        xyz_data = self._prepare_xyz_data(main_positions, side_positions)
        
        # Write .xyz file
        with open(file_path, 'w') as f:
            # First line: number of atoms
            f.write(f"{len(xyz_data)}\n")
            # Second line: comment
            f.write(f"{comment}\n")
            # Data lines: atom_type x y z
            for atom_type, x, y, z in xyz_data:
                f.write(f"{atom_type} {x:.6f} {y:.6f} {z:.6f}\n")
                
        logger.info(f"Saved .xyz file: {file_path}")
    
    def _prepare_xyz_data(self, main_positions: np.ndarray, side_positions: List[Optional[np.ndarray]]) -> List[Tuple[str, float, float, float]]:
        """Prepare data for .xyz file format."""
        xyz_data = []
        
        # Add main chain atoms
        for i, (amino_acid, position) in enumerate(zip(self.main_chain_sequence, main_positions)):
            xyz_data.append((amino_acid, position[0], position[1], position[2]))
            
        # Add side chain atoms
        if side_positions:
            for i, (amino_acid, position) in enumerate(zip(self.side_chain_sequence, side_positions)):
                if position is not None and amino_acid != EMPTY_SIDECHAIN_PLACEHOLDER:
                    xyz_data.append((amino_acid, position[0], position[1], position[2]))
                    
        return xyz_data
    
    def plot_3d_structure(self, 
                         main_positions: Optional[np.ndarray] = None,
                         side_positions: Optional[List[Optional[np.ndarray]]] = None,
                         title: str = "Protein 3D Structure",
                         show_grid: bool = True,
                         show_ticks: bool = True,
                         smooth_lines: bool = True,
                         save_path: Optional[str] = None) -> plt.Figure:
        """
        Create 3D plot of protein structure with smooth tetrahedral lattice.
        
        Args:
            main_positions: Main chain coordinates (uses cached if None)
            side_positions: Side chain coordinates (uses cached if None)
            title: Plot title
            show_grid: Whether to show grid
            show_ticks: Whether to show axis ticks
            smooth_lines: Whether to use smooth line interpolation
            save_path: Path to save plot (optional)
            
        Returns:
            matplotlib Figure object
        """
        if main_positions is None:
            main_positions = self._main_positions
        if side_positions is None:
            side_positions = self._side_positions
            
        if main_positions is None:
            raise ValueError("No coordinates available. Generate coordinates first.")
            
        # Create figure and 3D axis
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection='3d')
        
        # Plot main chain
        self._plot_main_chain(ax, main_positions, smooth_lines)
        
        # Plot side chains
        if side_positions:
            self._plot_side_chains(ax, main_positions, side_positions)
            
        # Customize plot
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_xlabel('X', fontsize=12)
        ax.set_ylabel('Y', fontsize=12)
        ax.set_zlabel('Z', fontsize=12)
        
        if not show_ticks:
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])
            
        ax.grid(show_grid)
        
        # Equal aspect ratio
        self._set_equal_aspect_3d(ax, main_positions, side_positions)
        
        # Add legend
        ax.legend()
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Saved plot: {save_path}")
            
        return fig
    
    def _plot_main_chain(self, ax, positions: np.ndarray, smooth_lines: bool):
        """Plot main chain with beads and connections."""
        x, y, z = positions[:, 0], positions[:, 1], positions[:, 2]
        
        # Plot connections
        if smooth_lines:
            ax.plot(x, y, z, 'b-', linewidth=2, alpha=0.7, label='Main Chain')
        else:
            for i in range(len(positions) - 1):
                ax.plot([x[i], x[i+1]], [y[i], y[i+1]], [z[i], z[i+1]], 
                       'b-', linewidth=2, alpha=0.7)
                
        # Plot beads
        colors = plt.cm.viridis(np.linspace(0, 1, len(positions)))
        scatter = ax.scatter(x, y, z, c=colors, s=100, alpha=0.8, edgecolors='black', linewidth=1)
        
        # Add amino acid labels
        for i, (pos, aa) in enumerate(zip(positions, self.main_chain_sequence)):
            ax.text(pos[0], pos[1], pos[2], f'{aa}{i}', fontsize=8)
    
    def _plot_side_chains(self, ax, main_positions: np.ndarray, side_positions: List[Optional[np.ndarray]]):
        """Plot side chains."""
        for i, (main_pos, side_pos, side_aa) in enumerate(zip(main_positions, side_positions, self.side_chain_sequence)):
            if side_pos is not None and side_aa != EMPTY_SIDECHAIN_PLACEHOLDER:
                # Connection from main chain to side chain
                ax.plot([main_pos[0], side_pos[0]], 
                       [main_pos[1], side_pos[1]], 
                       [main_pos[2], side_pos[2]], 
                       'r--', linewidth=1.5, alpha=0.6)
                
                # Side chain bead
                ax.scatter(side_pos[0], side_pos[1], side_pos[2], 
                          c='red', s=80, alpha=0.7, marker='s', 
                          edgecolors='darkred', linewidth=1, label='Side Chain' if i == 0 else "")
                
                # Side chain label
                ax.text(side_pos[0], side_pos[1], side_pos[2], f'{side_aa}', fontsize=7)
    
    def _set_equal_aspect_3d(self, ax, main_positions: np.ndarray, side_positions: Optional[List[Optional[np.ndarray]]]):
        """Set equal aspect ratio for 3D plot."""
        all_positions = [main_positions]
        
        if side_positions:
            valid_side_positions = [pos for pos in side_positions if pos is not None]
            if valid_side_positions:
                all_positions.extend(valid_side_positions)
                
        all_coords = np.vstack(all_positions)
        
        # Get the range of coordinates
        max_range = np.array([all_coords[:,0].max()-all_coords[:,0].min(),
                             all_coords[:,1].max()-all_coords[:,1].min(),
                             all_coords[:,2].max()-all_coords[:,2].min()]).max() / 2.0
        
        mid_x = (all_coords[:,0].max()+all_coords[:,0].min()) * 0.5
        mid_y = (all_coords[:,1].max()+all_coords[:,1].min()) * 0.5
        mid_z = (all_coords[:,2].max()+all_coords[:,2].min()) * 0.5
        
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)


def visualize_vqe_result(protein: Protein, 
                        vqe_result: str, 
                        unused_qubits: Optional[List[int]] = None,
                        output_dir: str = "",
                        filename: Optional[str] = None,
                        show_plot: bool = True,
                        save_plot: bool = False) -> ProteinFoldingVisualizer:
    """
    Convenience function to visualize VQE result in one call.
    
    Args:
        protein: Protein object
        vqe_result: VQE bitstring result
        unused_qubits: List of unused qubit indices
        output_dir: Directory for output files
        filename: Base filename for outputs
        show_plot: Whether to display the plot
        save_plot: Whether to save the plot
        
    Returns:
        ProteinFoldingVisualizer instance with computed coordinates
    """
    visualizer = ProteinFoldingVisualizer(protein)
    
    # Decode VQE result
    main_turns, side_turns = visualizer.decode_vqe_result(vqe_result, unused_qubits)
    
    # Generate coordinates
    main_pos, side_pos = visualizer.generate_coordinates(main_turns, side_turns)
    
    # Generate filename if not provided
    if filename is None:
        filename = f"protein_{''.join(protein.main_chain.sequence)}"
    
    # Save .xyz file
    try:
        visualizer.save_xyz_file(filename, path=output_dir, 
                               comment=f"VQE Result: {vqe_result}")
        logger.info(f"Generated .xyz file: {filename}.xyz")
    except Exception as e:
        logger.error(f"Failed to save .xyz file: {e}")
    
    # Create and optionally save plot
    try:
        fig = visualizer.plot_3d_structure(title=f"Protein Structure: {filename}")
        
        if save_plot:
            plot_path = os.path.join(output_dir, f"{filename}_3d.png")
            fig.savefig(plot_path, dpi=300, bbox_inches='tight')
            
        if show_plot:
            plt.show()
        else:
            plt.close(fig)
            
    except Exception as e:
        logger.error(f"Failed to create plot: {e}")
        
    return visualizer