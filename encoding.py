from enum import IntEnum, Enum
import numpy as np
import matplotlib.pyplot as plt


class SubLattice(IntEnum):
    A = 0
    B = 1


class Turn(Enum):
    DIR_0 = (1, 1, 1)
    DIR_1 = (1, -1, -1)
    DIR_2 = (-1, 1, -1)
    DIR_3 = (-1, -1, 1)


class Penalty(IntEnum):
    COLLISION = 100
    BACK = 1
    HYDROPHOBIC_HYDROPHOBIC = -1


class TetrahedralLattice:
    """
    A class for representing proteins on a tetrahedral lattice using FCC structure.
    
    Key consistency principle: All spatial measurements (bond lengths, move vectors, 
    and energy calculations) are scaled consistently with fcc_edge_length to ensure 
    that the relative geometry remains constant regardless of the absolute scale.
    
    The fundamental unit is the bond length, which equals:
    sqrt(0.25² + 0.25² + 0.25²) * fcc_edge_length = sqrt(3)/4 * fcc_edge_length
    """
    def __init__(self, fcc_edge_length=2.0, tolerance=1e-3):
        self.fcc_edge_length = fcc_edge_length
        self.tolerance = tolerance
        # Assuming the lattice will be limited to certain dimensions (we won't generate very long protein chains) - for now, we can
        # store all nodes, bonds and neighbors in memory
        self.nodes = []
        self.bonds = []
        self.cell_indices = []  # only for debugging purposes
        self.neighbors = {}

        self.move_vectors = None
        self._init_turn_vectors()

    def _get_bond_length(self):
        """
        Calculate the bond length for the FCC tetrahedral lattice.
        This is the distance between nearest neighbor positions in the base unit cell.
        """
        return np.linalg.norm([0.25, 0.25, 0.25]) * self.fcc_edge_length

    def _init_turn_vectors(self):
        raw = np.array(
            [turn.value for turn in Turn],
            dtype=float,
        )

        # Normalize the direction vectors
        normed = raw / np.linalg.norm(raw, axis=1)[:, None]
        
        # Use consistent bond length calculation
        bond_length = self._get_bond_length()
        
        # Scale the normalized direction vectors by the bond length
        self.move_vectors = normed * bond_length

    def generate_lattice(self, nx, ny, nz):
        # Base positions for the FCC lattice in a tetrahedral arrangement (8 nodes in one unit cell)
        # nx, ny, nz are the number of unit cells in each direction (NOT THE NUMBER OF NODES!!)

        base = (
            np.array(
                [
                    [0, 0, 0],
                    [0.25, 0.25, 0.25],
                    [0.5, 0.5, 0],
                    [0.75, 0.75, 0.25],
                    [0.5, 0, 0.5],
                    [0.75, 0.25, 0.75],
                    [0, 0.5, 0.5],
                    [0.25, 0.75, 0.75],
                ]
            )
            * self.fcc_edge_length
        )
        nodes = []
        cell_indices = []

        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    cell_index = (
                        ix * ny * nz + iy * nz + iz
                    )  # unique index for each cell
                    shift = np.array([ix, iy, iz]) * self.fcc_edge_length
                    for base_index, base_direction in enumerate(base):
                        pos = base_direction + shift
                        nodes.append(pos)

                        # debug purposes only
                        cell_indices.append(cell_index * len(base) + base_index)

        self.nodes = np.array(nodes)
        self.cell_indices = np.array(cell_indices)
        self._find_neighbors()

    def _find_neighbors(self):
        n = len(self.nodes)
        self.neighbors = {i: [] for i in range(n)}
        self.bonds = []
        
        # Use consistent bond length calculation
        bond_length = self._get_bond_length()
        
        for i in range(n):
            for j in range(i + 1, n):
                dist = np.linalg.norm(self.nodes[i] - self.nodes[j])
                if abs(dist - bond_length) < self.tolerance:
                    self.neighbors[i].append(j)
                    self.neighbors[j].append(i)
                    self.bonds.append((i, j))

    def _get_available_turns(self, sublattice):
        # TODO: Think over how to handle sublattice A and B
        # For now, let sublattice A be the original move vectors and sublattice B be the negative of them
        if sublattice == SubLattice.A:
            return self.move_vectors
        else:
            return -1 * self.move_vectors

    def generate_protein_path(self, beads, turn_sequence, starting_pos=None):
        if len(turn_sequence) != len(beads) - 1:
            raise ValueError("Turn sequence length must be equal to: len(beads) - 1")

        if starting_pos is None:
            starting_pos = np.array([0.0, 0.0, 0.0])
        positions = [starting_pos]

        for i, turn_idx in enumerate(turn_sequence):
            current_bead = beads[i]
            available_turns = self._get_available_turns(current_bead.sublattice)

            if turn_idx >= len(available_turns):
                raise ValueError(
                    f"Turn index {turn_idx} is out of range for {current_bead.sublattice} sublattice"
                )

            turning_move_vector = available_turns[turn_idx]
            new_position = positions[-1] + turning_move_vector
            positions.append(new_position)

        return np.array(positions)

    def compute_energy(self, positions, beads):
        """
        Compute the energy of a protein conformation based on:
        1. Hydrophobic-hydrophobic contacts (favorable, -1 energy)
        2. Collision detection (unfavorable, +100 energy)
        3. Backtracking detection (unfavorable, +1 energy)
        """
        energy = 0
        n = len(positions)
        
        # Use consistent bond length for distance calculations
        bond_length = self._get_bond_length()
        bond_length_squared = bond_length ** 2

        for i in range(n):
            for j in range(i + 1, n):
                if abs(i - j) > 1:
                    if beads[i].symbol == "H" and beads[j].symbol == "H":
                        # Check if H-H beads are at nearest neighbor distance
                        distance_squared = np.sum((positions[j] - positions[i]) ** 2)
                        if abs(distance_squared - bond_length_squared) < self.tolerance:
                            energy += Penalty.HYDROPHOBIC_HYDROPHOBIC
                            print("H-H found")

                    if np.allclose(positions[i], positions[j], atol=self.tolerance):
                        energy += Penalty.COLLISION
                        print(f"Collision detected between beads {i} and {j}")

        for i in range(2, n):
            if np.allclose(positions[i], positions[i - 2], atol=self.tolerance):
                energy += Penalty.BACK
                print(f"Backtracking detected at bead {i}")

        return energy

    def find_lowest_energy_conformation(
        self, beads, all_turn_sequences, starting_pos=[5.0, 5.0, 5.0]
    ):
        best_energy = float("inf")
        best_turns = None
        best_positions = None

        for turn_sequence in all_turn_sequences:
            try:
                positions = self.generate_protein_path(
                    beads, turn_sequence, starting_pos=np.array(starting_pos)
                )
                energy = self.compute_energy(positions, beads)

                if energy < best_energy:
                    best_energy = energy
                    best_turns = (
                        turn_sequence.copy()
                        if hasattr(turn_sequence, "copy")
                        else list(turn_sequence)
                    )
                    best_positions = positions.copy()
            except (ValueError, IndexError) as e:
                print(f"Skipping invalid conformation: {e}")
                continue

        return {
            "best_turns": best_turns,
            "best_energy": best_energy,
            "best_positions": best_positions,
        }

    def visualize_lattice(
        self,
        show_bonds=True,
        show_node_labels=True,
        protein_path=None,
        protein_sequence=None,
    ):
        fig = plt.figure(figsize=(14, 12))
        ax = fig.add_subplot(111, projection="3d")

        if show_bonds:
            for i, j in self.bonds:
                x = [self.nodes[i][0], self.nodes[j][0]]
                y = [self.nodes[i][1], self.nodes[j][1]]
                z = [self.nodes[i][2], self.nodes[j][2]]
                ax.plot(x, y, z, c="lightgray", alpha=0.4, linewidth=1, zorder=1)

        xs_all, ys_all, zs_all = self.nodes[:, 0], self.nodes[:, 1], self.nodes[:, 2]
        ax.scatter(
            xs_all,
            ys_all,
            zs_all,
            c="lightgray",
            s=60,
            alpha=0.4,
            label="Lattice nodes",
            zorder=2,
            depthshade=False,
        )

        if protein_path is not None:
            for k in range(len(protein_path) - 1):
                x = [protein_path[k][0], protein_path[k + 1][0]]
                y = [protein_path[k][1], protein_path[k + 1][1]]
                z = [protein_path[k][2], protein_path[k + 1][2]]
                ax.plot(x, y, z, c="red", linewidth=3, zorder=4)

            xs_p = protein_path[:, 0]
            ys_p = protein_path[:, 1]
            zs_p = protein_path[:, 2]
            ax.scatter(
                xs_p,
                ys_p,
                zs_p,
                c="green",
                s=200,
                edgecolors="black",
                label="Protein nodes",
                zorder=5,
            )

            if protein_sequence:
                for idx, (x, y, z, aa) in enumerate(
                    zip(xs_p, ys_p, zs_p, protein_sequence)
                ):
                    ax.text(
                        x,
                        y,
                        z + 0.2,
                        f"{aa}{idx}",
                        color="black",
                        fontsize=10,
                        ha="center",
                        zorder=6,
                    )

        if show_node_labels:
            for i, (x, y, z) in enumerate(self.nodes):
                ax.text(x, y, z, str(i), color="darkgray", fontsize=6, zorder=3)

        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title("Tetrahedral Lattice with Folded Protein Highlighted")
        if protein_path is not None:
            ax.legend()
        plt.tight_layout()
        plt.show()
        return fig

    def visualize_node_environment(self, node_id):
        """
        Visualize only one node and it's neighbors in 3D.
        """
        fig = plt.figure(figsize=(8, 7))
        ax = fig.add_subplot(111, projection="3d")
        x0, y0, z0 = self.nodes[node_id]
        ax.scatter([x0], [y0], [z0], c="red", s=150, label=f"Node {node_id}")
        colors = ["blue", "green", "orange", "purple"]
        for idx, neighbor in enumerate(self.neighbors[node_id]):
            xn, yn, zn = self.nodes[neighbor]
            ax.scatter(
                [xn], [yn], [zn], c=colors[idx % 4], s=100, label=f"Neighbor {neighbor}"
            )
            ax.plot([x0, xn], [y0, yn], [z0, zn], c=colors[idx % 4], linewidth=2)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title(f"Node {node_id} and its neighbors")
        ax.legend()
        plt.tight_layout()
        plt.show()
        return fig

    def print_neighbors(self):
        # IMPORTANT: Indexes of lattice nodes are not the same as bead indexes!!
        for i, neigh in self.neighbors.items():
            print(f"Node index: {i}: neighbors {neigh}")
