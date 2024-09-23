from triqs_maxent import *
import numpy as np
from h5 import HDFArchive
from triqs.gf import make_gf_from_fourier
from triqs.plot.mpl_interface import oplot, plt
from triqs_dft_tools.sumk_dft import *

h_field = 0.0
beta = 40
dft_filename = f"nsp"
use_blocks = True  # use bloc structure from DFT input
prec_mu = 0.0001  # precision of chemical potential
norb = 20

SK = SumkDFT(hdf_file=dft_filename + ".h5", use_dft_blocks=use_blocks, beta=beta)

with HDFArchive(dft_filename + ".h5") as ar:
    n_iterations = ar["DMFT_results"]["iteration_count"]
    G_latt_orb = ar["DMFT_results"]["Iterations"]["G_latt_orb_it" + str(n_iterations)]


results_latt_orb = {}
for name, giw in G_latt_orb:
    for i in range(norb):
        tm = TauMaxEnt()
        tm.set_G_iw(giw[i, i])
        tm.omega = HyperbolicOmegaMesh(omega_min=-20, omega_max=20, n_points=201)
        tm.alpha_mesh = LogAlphaMesh(alpha_min=1e-2, alpha_max=1e2, n_points=30)
        tm.set_error(1.0e-3)
        results_latt_orb[f"{name}_{i}"] = tm.run()


# HDF5 ファイルに MaxEnt 結果と omega を保存
with HDFArchive("results_maxent.h5", "w") as ar:
    for key, result in results_latt_orb.items():
        ar[key] = {
            "data": result.data,  # MaxEntResult のデータ部分
            "omega": result.omega,  # MaxEntResult の omega 部分
        }
