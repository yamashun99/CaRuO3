import numpy as np
import matplotlib.pyplot as plt
from extract_data import extract_data_from_h5


def plot_dmft_data(filename):
    data = extract_data_from_h5(filename)

    n_iterations = data["n_iterations"]
    iterations = range(1, n_iterations + 1)

    # プロットの準備 (sharex=Trueを追加して横軸を共有)
    fig, axs = plt.subplots(5, 1, figsize=(8, 15), sharex=True)

    # 化学ポテンシャルのプロット
    axs[0].plot(iterations, data["chemical_potentials"], marker="o")
    axs[0].set_ylabel("Chemical Potential")
    axs[0].set_title("Chemical Potential over iterations")
    axs[0].grid(True)

    # GimpとGlocの差のプロット
    axs[1].plot(iterations, data["G_diff"], marker="o")
    axs[1].set_ylabel("Diff between Gloc and Gimp")
    axs[1].set_title("Difference between Gloc and Gimp over iterations")
    axs[1].grid(True)
    axs[1].set_yscale("log")

    # Gimpの差分のプロット
    axs[2].plot(iterations, data["Gimp_diff"], marker="o")
    axs[2].set_ylabel("Diff between Gimp(n) and Gimp(n-1)")
    axs[2].set_title("Difference between Gimp(n) and Gimp(n-1) over iterations")
    axs[2].grid(True)
    axs[2].set_yscale("log")

    # 密度行列のプロット
    axs[3].plot(iterations, data["density_matrix_values"], marker="o")
    axs[3].set_ylabel("Density Matrix Value")
    axs[3].set_title("Density Matrix over iterations")
    axs[3].grid(True)

    # 軌道ごとの密度行列のプロット
    real_density_data = {}
    for data_per_iter in data["orbital_resolved_dms"]:
        for key in data_per_iter:
            if key not in real_density_data:
                real_density_data[key] = []
            real_density_data[key].append(np.real(data_per_iter[key]))

    for key, values in real_density_data.items():
        axs[4].plot(iterations, values, label=key, marker="o")

    axs[4].set_xlabel("Iteration")
    axs[4].set_ylabel("Orbital Resolved Density (Real Part)")
    axs[4].set_title("Orbital Resolved Density Matrix over Iterations")
    axs[4].grid(True)
    axs[4].legend()

    plt.tight_layout()
    plt.show()
    plt.savefig("convergence.pdf")


plot_dmft_data("nsp.h5")
