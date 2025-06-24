from encoding import TetrahedralLattice, all_turn_combinations
import matplotlib.pyplot as plt


def main() -> None:
    main_chain = ["H", "P", "P", "H", "P", "H"]

    lattice = TetrahedralLattice(30, 30, 30)
    all_turns = all_turn_combinations(len(main_chain))

    for turns in all_turns:
        try:
            qubit_string = lattice.encode_turn_sequence(turns)

            positions = lattice.generate_positions((10, 10, 10), turns)
            xs = [p[0] for p in positions]
            ys = [p[1] for p in positions]
            zs = [p[2] for p in positions]

            with open("encoded_turns.txt", "a+") as f:
                f.write(f"x:{xs} y:{ys} z:{zs} | string: {qubit_string}\n")

            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")
            ax.plot(xs, ys, zs, marker="o")
            ax.set_title("HP conformations on tetrahedral lattice")
            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            ax.set_zlabel("Z")

            plt.show()
            # plt.savefig("plots/plot_" + qubit_string + ".png")
            # plt.close(fig)
        except ValueError as e:
            print(f"Error: {e}")


if __name__ == "__main__":
    main()
