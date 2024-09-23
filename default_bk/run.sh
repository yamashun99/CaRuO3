#!/bin/bash

# 配列の定義
x_values=(5 7 10 20)
y_values=(0.02 0.01 0.005)
z_values=(0.25 0.5)

# ループ処理
for z in "${z_values[@]}"; do
    for x in "${x_values[@]}"; do
        for y in "${y_values[@]}"; do
            # ディレクトリ作成
            dir_name="b${x}-h${y}-U3.0-J${z}"
            mkdir "$dir_name"

            # ファイルコピー
            cp ./nsp.h5 "$dir_name"
            cp ./extract_data.py "$dir_name"
            cp ./dmft.py "$dir_name"

            # ディレクトリに移動
            cd "$dir_name"

            # dmft.pyの内容を書き換える
            sed -i "s/h_field = .*/h_field = ${y}/" dmft.py
            sed -i "s/beta = .*/beta = ${x}/" dmft.py
            sed -i "s/J = .*/J = ${z}/" dmft.py

            # スクリプトを実行
            time mpirun python dmft.py | tee dmft.out

            # 元のディレクトリに戻る
            cd ..
        done
    done
done
