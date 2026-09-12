#!/usr/bin/env bash
#
# pick_every_n.sh
#
# 入力ディレクトリ内のファイルをファイル名順(=時刻順のはず)に並べ、
# N枚に1枚だけを別ディレクトリへコピーする。
# ファイル名にスペースが入っていても安全に扱える(find -print0 + 配列)。
#
# 使い方:
#   ./pick_every_n.sh 入力ディレクトリ 出力ディレクトリ N [glob パターン]
#
# 例:
#   # "Light_rz cas_*.fit*" にマッチするファイルから20枚に1枚を間引く
#   ./pick_every_n.sh ./rz_casG ./rz_casG_sparse 20 "Light_rz cas_*.fit*"
#
#   # パターン省略時は全ファイルが対象
#   ./pick_every_n.sh ./rz_casG ./rz_casG_sparse 20

set -euo pipefail

if [ $# -lt 3 ]; then
    echo "使い方: $0 入力ディレクトリ 出力ディレクトリ N [globパターン]" >&2
    exit 1
fi

input_dir=$1
output_dir=$2
n=$3
pattern=${4:-*}

if [ ! -d "$input_dir" ]; then
    echo "エラー: 入力ディレクトリが見つかりません: $input_dir" >&2
    exit 1
fi

if ! [[ "$n" =~ ^[0-9]+$ ]] || [ "$n" -lt 1 ]; then
    echo "エラー: N は1以上の整数で指定してください(指定値: $n)" >&2
    exit 1
fi

mkdir -p "$output_dir"

# ファイル名にスペースが入っていても壊れないよう、NUL区切りで受け取る。
# find の結果はディレクトリ内をそのまま列挙するので、ファイル名順に
# ソートし直す(時刻がファイル名に埋め込まれている前提)。
mapfile -d '' -t files < <(find "$input_dir" -maxdepth 1 -type f -name "$pattern" -print0 | sort -z)

total=${#files[@]}
if [ "$total" -eq 0 ]; then
    echo "エラー: '$pattern' に一致するファイルが $input_dir にありません。" >&2
    exit 1
fi

copied=0
for ((i = 0; i < total; i += n)); do
    src=${files[$i]}
    cp -- "$src" "$output_dir/"
    copied=$((copied + 1))
done

echo "完了: ${total}枚中 ${copied}枚を ${n}枚に1枚の割合でコピーしました → $output_dir"
