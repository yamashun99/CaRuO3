from itertools import *
import numpy as np
import triqs.utility.mpi as mpi
from h5 import *
from triqs.gf import *
from triqs_dft_tools.sumk_dft import *
from triqs_dft_tools.sumk_dft_tools import *
from triqs.operators.util.hamiltonians import *
from triqs.operators.util.U_matrix import *
from triqs_cthyb import *
import warnings
from triqs.plot.mpl_interface import oplot, plt

warnings.filterwarnings("ignore", category=FutureWarning)
filename = f"nsp"
with HDFArchive(filename + ".h5", "r") as ar:
    dmft_results = ar["DMFT_results"]["Iterations"]
    n_iterations = ar["DMFT_results"]["iteration_count"]
    Gimp = dmft_results["Gimp_it" + str(n_iterations - 1)]
SK = SumkDFTTools(hdf_file=filename + ".h5", use_dft_blocks=False, mesh=Gimp.mesh)

# We analyze the block structure of the Hamiltonian
Sigma = SK.block_structure.create_gf(mesh=Gimp.mesh)

SK.put_Sigma([Sigma])

# Setup CTQMC Solver
n_orb = SK.corr_shells[0]["dim"]
spin_names = ["up", "down"]
# Print some information on the master node
Sigma_iw = 0
with HDFArchive(filename + ".h5", "r") as ar:
    print(ar)
    SK.block_structure = ar["DMFT_input"]["sumk_block_structure"]
    iteration_offset = ar["DMFT_results"]["iteration_count"] + 1
    print(("offset", iteration_offset))
    Sigma_iw = ar["DMFT_results"]["Iterations"]["Sigma_it" + str(iteration_offset - 1)]
    SK.dc_imp = ar["DMFT_results"]["Iterations"]["dc_imp" + str(iteration_offset - 1)]
    SK.dc_energ = ar["DMFT_results"]["Iterations"][
        "dc_energ" + str(iteration_offset - 1)
    ]
    SK.chemical_potential = ar["DMFT_results"]["Iterations"][
        "chemical_potential" + str(iteration_offset - 1)
    ].real
SK.put_Sigma(Sigma_imp=[Sigma_iw])

ikarray = numpy.array(list(range(SK.n_k)))

gf_csc = Gf(mesh=SK.mesh, target_shape=(SK.lattice_gf(0)["up"].target_shape))
G_latt_orb = BlockGf(
    name_list=["up", "down"], block_list=[gf_csc, gf_csc], make_copies=True
)


for ik in mpi.slice_array(ikarray):
    G_latt_KS = SK.lattice_gf(ik=ik) * SK.bz_weights[ik]
    for bname, gf in G_latt_orb:
        gf += G_latt_KS[bname]
nup = np.sum(np.diag(G_latt_orb.density()["up"]))
print(f"nup:{nup}")
ndown = np.sum(np.diag(G_latt_orb.density()["down"]))
print(f"ndown:{ndown}")

with HDFArchive(filename + ".h5", "a") as ar:
    ar["DMFT_results"]["Iterations"][
        "G_latt_orb_it" + str(iteration_offset - 1)
    ] = G_latt_orb
