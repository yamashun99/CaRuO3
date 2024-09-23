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

            # ファイルコピー
            cp ./local_lattice.py "$dir_name"

            # ディレクトリに移動
            cd "$dir_name"

            # スクリプトを実行
            python local_lattice.py

            # 元のディレクトリに戻る
            cd ..
        done
    done
done
