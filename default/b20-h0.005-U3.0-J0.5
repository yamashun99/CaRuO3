from itertools import *
import numpy as np
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

ar = HDFArchive(filename + ".h5", "a")
# HDF5ファイルを最初に読み込む（読み取り専用）
dmft_results = ar["DMFT_results"]["Iterations"]
n_iterations = ar["DMFT_results"]["iteration_count"]
Gimp = dmft_results["Gimp_it" + str(n_iterations - 1)]

# SumkDFTToolsの初期化
SK = SumkDFTTools(hdf_file=filename + ".h5", use_dft_blocks=False, mesh=Gimp.mesh)

# ブロック構造の自己エネルギーの計算
Sigma = SK.block_structure.create_gf(mesh=Gimp.mesh)
SK.put_Sigma([Sigma])

# HDF5ファイルの操作（書き込みモードで開く）
Sigma_iw = 0
block_structure = None
iteration_offset = None
block_structure = ar["DMFT_input"]["sumk_block_structure"]
iteration_offset = ar["DMFT_results"]["iteration_count"] + 1
print(("offset", iteration_offset))
Sigma_iw = ar["DMFT_results"]["Iterations"]["Sigma_it" + str(iteration_offset - 1)]
SK.dc_imp = ar["DMFT_results"]["Iterations"]["dc_imp" + str(iteration_offset - 1)]
SK.dc_energ = ar["DMFT_results"]["Iterations"]["dc_energ" + str(iteration_offset - 1)]
SK.chemical_potential = ar["DMFT_results"]["Iterations"][
    "chemical_potential" + str(iteration_offset - 1)
].real


# ブロック構造の確認と設定
if block_structure:
    SK.block_structure = block_structure
else:
    G = SK.extract_G_loc(transform_to_solver_blocks=False)
    SK.analyse_block_structure_from_gf(G, threshold=1e-3)

# 自己エネルギーをセット
SK.put_Sigma(Sigma_imp=[Sigma_iw])

ikarray = numpy.array(list(range(SK.n_k)))
# 格子グリーン関数の計算
gf_csc = Gf(mesh=SK.mesh, target_shape=(SK.lattice_gf(0)["up"].target_shape))
G_latt_orb = BlockGf(
    name_list=["up", "down"], block_list=[gf_csc, gf_csc], make_copies=True
)

# k点ごとの格子グリーン関数の計算
for ik in ikarray:
    G_latt_KS = SK.lattice_gf(ik=ik) * SK.bz_weights[ik]
    for bname, gf in G_latt_orb:
        gf += G_latt_KS[bname]


## マスターノードで結果をHDF5ファイルに書き込む
ar["DMFT_results"]["Iterations"][
    "G_latt_orb_it" + str(iteration_offset - 1)
] = G_latt_orb
