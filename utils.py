import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def draw_chain(positions, main_chain, lattice=None):  

    xs = [p[0] for p in positions]
    ys = [p[1] for p in positions]
    zs = [p[2] for p in positions]
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    
    # Only draw tetrahedral lattice connections if lattice is provided
    if lattice is not None:
        # Get all unique positions in our protein
        unique_positions = list(set(tuple(p) for p in positions))
        
        print("=== Ilość sąsiadów dla każdej pozycji ===")
        
        # For each position, draw lines only to neighbors that are also in the protein
        drawn_connections = set()  # To avoid drawing the same connection twice
        
        for i, pos in enumerate(unique_positions):
            neighbors = lattice.get_neighbors(pos)
            protein_neighbors = [neighbor for neighbor in neighbors if neighbor in unique_positions]
            
            print(f"Pozycja {i}: {pos}")
            print(f"  Wszyscy sąsiedzi tetraedryczni: {len(neighbors)}")
            print(f"  Sąsiedzi w białku: {len(protein_neighbors)}")
            print(f"  Lista sąsiadów: {protein_neighbors}")
            
            for neighbor in neighbors:
                # Only draw connection if the neighbor is also a position in our protein
                if neighbor in unique_positions:
                    # Create a sorted tuple to avoid duplicate connections
                    connection = tuple(sorted([pos, neighbor]))
                    if connection not in drawn_connections:
                        ax.plot([pos[0], neighbor[0]], 
                               [pos[1], neighbor[1]], 
                               [pos[2], neighbor[2]], 
                               linestyle="--", color="lightblue", alpha=1, linewidth=1.5)
                        drawn_connections.add(connection)
        print()
    
    # Draw the actual chain connections with solid lines
    ax.plot(xs, ys, zs, marker="o", markersize=12, linestyle="-", 
           color="darkgray", linewidth=2)

    # Annotate each point with "H" or "P"
    for x, y, z, label in zip(xs, ys, zs, main_chain):
        ax.text(x, y, z, label, fontsize=14, ha='center', va='center', color='red' if label == 'H' else 'blue')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.grid(True)

    ax.set_title("Tetrahedral Lattice")
    plt.show()