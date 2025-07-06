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
    def __init__(self, fcc_edge_length=2.0, tolerance=1e-3):
        self.fcc_edge_length = fcc_edge_length
        self.tolerance = tolerance
        # Assuming the lattice will be limited to certain dimensions (we won't generate very long protein chains) - for now, we can
        # store all nodes, bonds and neighbors in memory
        self.nodes = []
        self.bonds = []
        self.neighbors = {}

    def generate_lattice(self, nx, ny, nz):
        # Base positions for the FCC lattice in a tetrahedral arrangement (8 nodes in one unit cell)
        # nx, ny, nz are the number of unit cells in each direction (NOT THE NUMBER OF NODES!!)

        base = np.array([
            [0, 0, 0],
            [0.25, 0.25, 0.25],
            [0.5, 0.5, 0],
            [0.75, 0.75, 0.25],
            [0.5, 0, 0.5],
            [0.75, 0.25, 0.75],
            [0, 0.5, 0.5],
            [0.25, 0.75, 0.75]
        ])
        nodes = []
        
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    shift = np.array([ix, iy, iz])
                    for b in base:
                        pos = (b + shift) * self.fcc_edge_length
                        nodes.append(pos)
        
        self.nodes = np.array(nodes)
        self._find_neighbors()

    def _find_neighbors(self):
        n = len(self.nodes)
        self.neighbors = {i: [] for i in range(n)}
        self.bonds = []
        # TODO
        # Check the bond length for tetrahedral lattice
        # For FCC, the bond length should be: sqrt(3) * edge_length / 2 (Pythagorean theorem)

        bond_length = np.linalg.norm(np.array([0.25, 0.25, 0.25]) * self.fcc_edge_length)
        for i in range(n):
            for j in range(i+1, n):
                dist = np.linalg.norm(self.nodes[i] - self.nodes[j])
                if abs(dist - bond_length) < self.tolerance:
                    self.neighbors[i].append(j)
                    self.neighbors[j].append(i)
                    self.bonds.append((i, j))

    def _get_available_turns(self, sublattice):
        if sublattice == SubLattice.A:
            return [np.array(turn.value) for turn in Turn]
        else:
            return [-np.array(turn.value) for turn in Turn]

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
                raise ValueError(f"Turn index {turn_idx} is out of range for {current_bead.sublattice} sublattice")
            
            turning_move_vector = available_turns[turn_idx]
            new_position = positions[-1] + turning_move_vector
            positions.append(new_position)
        
        return np.array(positions)

    def visualize_lattice(self, show_bonds=True, show_node_labels=True, protein_path=None, protein_sequence=None):
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        
        for i, (x, y, z) in enumerate(self.nodes):
            ax.scatter([x], [y], [z], c='gray', s=80, alpha=0.7)

        if show_node_labels:
            for i, (x, y, z) in enumerate(self.nodes):
                ax.text(x, y, z, str(i), color='black', fontsize=6)
        
        if show_bonds:
            for i, j in self.bonds:
                x = [self.nodes[i][0], self.nodes[j][0]]
                y = [self.nodes[i][1], self.nodes[j][1]]
                z = [self.nodes[i][2], self.nodes[j][2]]
                ax.plot(x, y, z, c='gray', alpha=0.3, linewidth=1)
        
        if protein_path is not None:
            xs = protein_path[:, 0]
            ys = protein_path[:, 1]
            zs = protein_path[:, 2]
            ax.plot(xs, ys, zs, 'ro-', linewidth=2, markersize=5, label='Folded protein sequence')
            
            if protein_sequence:
                for i, (x, y, z, aa) in enumerate(zip(xs, ys, zs, protein_sequence)):
                    color = 'red' if aa == 'H' else 'blue'
                    ax.text(x, y, z + 0.2, f'{aa}{i}', color=color, fontsize=12, 
                           ha='center', va='center', weight='bold')
        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Tetrahedral Lattice Visualization')
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
        ax = fig.add_subplot(111, projection='3d')
        x0, y0, z0 = self.nodes[node_id]
        ax.scatter([x0], [y0], [z0], c='red', s=150, label=f'Node {node_id}')
        colors = ['blue', 'green', 'orange', 'purple']
        for idx, neighbor in enumerate(self.neighbors[node_id]):
            xn, yn, zn = self.nodes[neighbor]
            ax.scatter([xn], [yn], [zn], c=colors[idx%4], s=100, label=f'Neighbor {neighbor}')
            ax.plot([x0, xn], [y0, yn], [z0, zn], c=colors[idx%4], linewidth=2)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(f'Node {node_id} and its neighbors')
        ax.legend()
        plt.tight_layout()
        plt.show()
        return fig

    def print_neighbors(self):
        # IMPORTANT: Indexes of lattice nodes are not the same as bead indexes!!
        for i, neigh in self.neighbors.items():
            print(f"Node index: {i}: neighbors {neigh}")

