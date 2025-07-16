from typing import Type, Callable, Any
from encoding import TetrahedralLattice
from numpy.typing import NDArray
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

class Visualizer:
    def __init__(self) -> None:
        self._registry = {}
        self.register(TetrahedralLattice, self._visualize_lattice)

    def register(self, obj_type: Type, func: Callable) -> None:
        self._registry[obj_type] = func

    def visualize(self, obj: Any, **kwargs) -> None:
        obj_type = type(obj)
        if obj_type in self._registry:
            self._registry[obj_type](obj, **kwargs)
        else:
            raise TypeError(f"No visualization function registered for {obj_type.__name__}")

    def _visualize_lattice(
        self,
        obj,
        show_bonds: bool = True,
        show_node_labels: bool = True,
        protein_path: NDArray[np.float64] | None = None,
        protein_sequence: list[str] | None = None,
    ) -> Figure:
        fig = plt.figure(figsize=(14, 12))
        ax = fig.add_subplot(111, projection="3d")

        if show_bonds:
            for i, j in obj.bonds:
                x = [obj.nodes[i][0], obj.nodes[j][0]]
                y = [obj.nodes[i][1], obj.nodes[j][1]]
                z = [obj.nodes[i][2], obj.nodes[j][2]]
                ax.plot(x, y, z, c="lightgray", alpha=0.4, linewidth=1, zorder=1)

        xs_all, ys_all, zs_all = obj.nodes[:, 0], obj.nodes[:, 1], obj.nodes[:, 2]
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
            for i, (x, y, z) in enumerate(obj.nodes):
                ax.text(x, y, z, str(i), color="darkgray", fontsize=6, zorder=3)

        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")  # type: ignore
        ax.set_title("Tetrahedral Lattice with Folded Protein Highlighted")
        if protein_path is not None:
            ax.legend()
        plt.tight_layout()
        plt.show()
        return fig