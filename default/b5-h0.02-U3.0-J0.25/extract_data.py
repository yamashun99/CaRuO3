import numpy as np
from h5 import HDFArchive
from triqs.gf import *


def extract_data_from_h5(filename):
    with HDFArchive(filename, "r") as ar:
        dmft_results = ar["DMFT_results"]
        n_iterations = ar["DMFT_results"]["iteration_count"]

        # 化学ポテンシャルの取得
        chemical_potentials = []
        for it in range(n_iterations):
            iteration_key = f"chemical_potential{it}"
            if iteration_key in dmft_results["Iterations"]:
                chemical_potentials.append(dmft_results["Iterations"][iteration_key])

        # 化学ポテンシャルの差の取得
        mu_diff = []
        for it in range(1, n_iterations + 1):
            iteration_key = f"chemical_potential{it}"
            iteration_key_prev = f"chemical_potential{it-1}"
            if iteration_key in dmft_results["Iterations"]:
                mu_diff.append(
                    dmft_results["Iterations"][iteration_key]
                    - dmft_results["Iterations"][iteration_key_prev]
                )

        # GimpとGlocの差のプロット
        G_diff = []
        for it in range(1, n_iterations + 1):
            Gimp_it = ar["DMFT_results"]["Iterations"]["Gimp_it" + str(it)]
            Gloc_it = ar["DMFT_results"]["Iterations"]["Gloc_it" + str(it)]

            diff = max_G_diff(Gimp_it, Gloc_it)
            G_diff.append(diff)

        # Gimpの差分計算
        Gimp_diff = []
        for it in range(1, n_iterations + 1):
            Gimp_it = dmft_results["Iterations"]["Gimp_it" + str(it)]
            Gimp_prev = dmft_results["Iterations"]["Gimp_it" + str(it - 1)]
            diff = max_G_diff(Gimp_it, Gimp_prev)
            Gimp_diff.append(diff)

        # 密度行列の計算
        density_matrix_values = []
        for it in range(n_iterations):
            density_matrix_it = ar["DMFT_results"]["Iterations"][
                "Gloc_it" + str(it)
            ].density()
            total_real_sum = sum(
                float(np.real(np.sum(density_matrix_it[key])))
                for key in density_matrix_it
            )
            density_matrix_values.append(total_real_sum)

        # 密度行列の差の計算
        density_matrix_diff = []
        for it in range(1, n_iterations + 1):
            density_matrix_it = ar["DMFT_results"]["Iterations"][
                "Gloc_it" + str(it)
            ].density()
            density_matrix_prev = ar["DMFT_results"]["Iterations"][
                "Gloc_it" + str(it - 1)
            ].density()
            total_real_sum_it = sum(
                float(np.real(np.sum(density_matrix_it[key])))
                for key in density_matrix_it
            )
            total_real_sum_prev = sum(
                float(np.real(np.sum(density_matrix_prev[key])))
                for key in density_matrix_prev
            )
            density_matrix_diff.append(total_real_sum_it - total_real_sum_prev)

        # 軌道ごとの密度行列の取得
        orbital_resolved_dms = []
        for it in range(n_iterations):
            orbital_resolved_dm_it = ar["DMFT_results"]["Iterations"][
                "density_matrix_it" + str(it)
            ]
            flattened_array = {}
            for key in orbital_resolved_dm_it:
                flattened_array[key] = np.real(np.diag(orbital_resolved_dm_it[key]))
            orbital_resolved_dms.append(flattened_array)

        # 軌道ごとの密度行列の差の計算
        orbital_resolved_dm_diff = []
        for it in range(1, n_iterations + 1):
            orbital_resolved_dm_it = ar["DMFT_results"]["Iterations"][
                "density_matrix_it" + str(it)
            ]
            orbital_resolved_dm_prev = ar["DMFT_results"]["Iterations"][
                "density_matrix_it" + str(it - 1)
            ]
            diff = []
            for key in orbital_resolved_dm_it:
                diff.append(
                    (
                        np.diag(
                            np.real(orbital_resolved_dm_it[key])
                            - np.real(orbital_resolved_dm_prev[key])
                        )
                    )
                )
            orbital_resolved_dm_diff.append(diff)

    return {
        "chemical_potentials": chemical_potentials,
        "G_diff": G_diff,
        "Gimp_diff": Gimp_diff,
        "density_matrix_values": density_matrix_values,
        "orbital_resolved_dms": orbital_resolved_dms,
        "n_iterations": n_iterations,
        "mu_diff": mu_diff,
        "density_matrix_diff": density_matrix_diff,
        "orbital_resolved_dm_diff": orbital_resolved_dm_diff,
    }


def max_G_diff(G1, G2, norm_temp=True):
    if isinstance(G1, BlockGf):
        diff = 0.0
        for block, gf in G1:
            diff += max_G_diff(G1[block], G2[block], norm_temp)
        return diff

    assert G1.mesh == G2.mesh, "mesh of two input Gfs does not match"
    assert G1.target_shape == G2.target_shape, "can only compare Gfs with same shape"

    if type(G1.mesh) is MeshImFreq:
        offset = np.diag(np.diag(G1.data[-1, :, :].real - G2.data[-1, :, :].real))
    else:
        offset = 0.0

    norm_grid = abs(
        np.linalg.norm(G1.data - G2.data - offset, axis=tuple(range(1, G1.data.ndim)))
    )
    norm = np.linalg.norm(norm_grid, axis=0)

    if type(G1.mesh) is MeshImFreq:
        norm = np.linalg.norm(norm_grid, axis=0) / np.sqrt(G1.mesh.beta)
    elif type(G1.mesh) is MeshImTime:
        norm = np.linalg.norm(norm_grid, axis=0) * np.sqrt(G1.mesh.beta / len(G1.mesh))
    elif type(G1.mesh) is MeshReFreq:
        norm = np.linalg.norm(norm_grid, axis=0) / np.sqrt(len(G1.mesh))
    else:
        raise ValueError("MeshReTime is not implemented")

    if type(G1.mesh) in (MeshImFreq, MeshImTime) and norm_temp:
        norm = norm / np.sqrt(G1.mesh.beta)

    return norm
