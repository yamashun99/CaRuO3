import os
import numpy as np
from h5 import HDFArchive
from triqs.gf import make_gf_from_fourier
import json

# ループで使用する値のリスト
x_values = [5, 7, 10, 20]
y_values = [0.02, 0.01, 0.005]
z_values = [0.25, 0.5]

# ループ処理
for z in z_values:
    for x in x_values:
        for y in y_values:
            # ディレクトリ名の作成
            dir_name = f"b{x}-h{y}-U3.0-J{z}"

            # dft_filenameの作成
            dft_filename = f"{dir_name}/nsp"

            # HDFファイルの読み込み
            try:
                with HDFArchive(dft_filename + ".h5") as ar:
                    n_iterations = ar["DMFT_results"]["iteration_count"]
                    G_latt_orb = ar["DMFT_results"]["Iterations"][
                        "G_latt_orb_it" + str(n_iterations)
                    ]

                    # n_upとn_downの計算
                    n_up = np.sum(np.real(np.diag(G_latt_orb["up"].density())))
                    n_down = np.sum(np.real(np.diag(G_latt_orb["down"].density())))

                    # 結果をJSONファイルに保存
                    result_data = {"n_up": n_up, "n_down": n_down}
                    json_file_path = os.path.join(dir_name, "n_up_down.json")
                    with open(json_file_path, "w") as json_file:
                        json.dump(result_data, json_file, indent=4)
                    print(f"Results saved in {json_file_path}")

            except Exception as e:
                print(f"Error processing {dft_filename}: {e}")
