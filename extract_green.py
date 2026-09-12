#!/usr/bin/env python3
"""
extract_green.py

Seestar等のOSCセンサーで撮影・スタックされたFITSから、
PixInsightのデベイヤー(補間)+XISF変換を経由せずに、直接
G(緑)チャンネルだけを抽出する。

入力データの形が2種類あり得るので、自動判別する:

1. 生のBayerモノクロ画像 (2次元, 例: 1920x1080)
   -> 2x2 Bayerタイルの2つのG画素を平均する「Superpixel法」を使う。
      補間を使わない実測値ベースなのでSNRも良く測光向き。
      出力は縦横半分の解像度になる。

2. すでにデベイヤー済みのRGBキューブ (3次元, 例: 3x1920x1080)
   -> Seestarの「Stacked_...fit」はこちらのケースが多い
      (Seestar自身が内部でデベイヤーしてRGB 3プレーンで保存する)。
      この場合は単純にG面(通常インデックス1)を取り出すだけでよく、
      解像度はそのまま維持される。

使い方:
    python3 extract_green.py input.fits output_G.fits [--pattern RGGB]
    python3 extract_green.py input.fits output_G.fits --show-shape

--show-shape を付けると、実際のデータ形状・軸の並びだけを表示して終了する
(手元のファイルがどちらのケースか確認したいときに使う)。

--pattern は2次元(生Bayer)の場合のみ必要。省略時はヘッダーの
BAYERPAT キーワードを自動で読み取る。
"""

import argparse
import sys

import numpy as np
from astropy.io import fits


def superpixel_green(data: np.ndarray, pattern: str) -> np.ndarray:
    """2x2 Bayerタイルの2つのG画素を平均し、G画像(半解像度)を返す。"""
    pattern = pattern.upper()

    layouts = {
        "RGGB": {(0, 0): "R", (0, 1): "G", (1, 0): "G", (1, 1): "B"},
        "BGGR": {(0, 0): "B", (0, 1): "G", (1, 0): "G", (1, 1): "R"},
        "GRBG": {(0, 0): "G", (0, 1): "R", (1, 0): "B", (1, 1): "G"},
        "GBRG": {(0, 0): "G", (0, 1): "B", (1, 0): "R", (1, 1): "G"},
    }

    if pattern not in layouts:
        raise ValueError(
            f"未対応のBayerパターン: {pattern} "
            f"(対応: {', '.join(layouts.keys())})"
        )

    layout = layouts[pattern]
    g_positions = [pos for pos, c in layout.items() if c == "G"]

    h, w = data.shape
    h2, w2 = h // 2, w // 2

    acc = np.zeros((h2, w2), dtype=np.float64)
    for dr, dc in g_positions:
        acc += data[dr : dr + 2 * h2 : 2, dc : dc + 2 * w2 : 2].astype(np.float64)
    acc /= len(g_positions)

    if np.issubdtype(data.dtype, np.integer):
        return np.round(acc).astype(data.dtype)
    return acc.astype(data.dtype)


def extract_from_rgb_cube(data: np.ndarray) -> np.ndarray:
    """すでにR,G,B 3プレーンになっているキューブからG面を取り出す。"""
    # サイズが3(またはRGBAで4)の軸をカラー軸とみなす
    color_axis = None
    for axis, size in enumerate(data.shape):
        if size in (3, 4):
            color_axis = axis
            break

    if color_axis is None:
        raise ValueError(
            f"3次元データですが、サイズ3または4の軸(カラーチャンネル)が"
            f"見つかりません: shape={data.shape}"
        )

    # R,G,B(,A)の順を仮定し、インデックス1(G)を取り出す
    green = np.take(data, indices=1, axis=color_axis)
    return green


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="入力FITS")
    parser.add_argument("output", nargs="?", help="出力FITS (Gチャンネルのみ)")
    parser.add_argument(
        "--pattern",
        help="Bayerパターンを明示指定 (RGGB/BGGR/GRBG/GBRG)。"
        "2次元の生Bayerデータの場合のみ使用。省略時はヘッダーのBAYERPATを使う。",
    )
    parser.add_argument(
        "--show-shape",
        action="store_true",
        help="データ形状(shape)とヘッダーの主要キーワードだけ表示して終了する",
    )
    args = parser.parse_args()

    with fits.open(args.input, memmap=False) as hdul:
        hdu = hdul[0]
        data = hdu.data.copy() if hdu.data is not None else None
        header = hdu.header.copy()

    if data is None:
        print("エラー: プライマリHDUに画像データがありません。", file=sys.stderr)
        return 1

    if args.show_shape:
        print(f"shape: {data.shape}  dtype: {data.dtype}")
        for key in ("NAXIS", "NAXIS1", "NAXIS2", "NAXIS3", "BAYERPAT", "COLORTYP"):
            if key in header:
                print(f"{key} = {header[key]}")
        return 0

    if args.output is None:
        print("エラー: 出力ファイル名を指定してください。", file=sys.stderr)
        return 1

    if data.ndim == 2:
        pattern = args.pattern or header.get("BAYERPAT") or header.get("COLORTYP")
        if not pattern:
            print(
                "エラー: 2次元(生Bayer)データですが、Bayerパターンが"
                "指定されておらず、ヘッダーからも判定できませんでした。"
                "--pattern RGGB のように明示してください。",
                file=sys.stderr,
            )
            return 1

        green = superpixel_green(data, str(pattern))

        for key in ("NAXIS1", "NAXIS2"):
            if key in header:
                del header[key]
        for key in ("CRPIX1", "CRPIX2"):
            if key in header:
                header[key] = header[key] / 2.0
        for key in ("CDELT1", "CDELT2"):
            if key in header:
                header[key] = header[key] * 2.0
        for key in ("CD1_1", "CD1_2", "CD2_1", "CD2_2"):
            if key in header:
                header[key] = header[key] * 2.0

        header["HISTORY"] = f"Green channel extracted via superpixel ({pattern})"
        header["BAYERPAT"] = "NONE"

    elif data.ndim == 3:
        green = extract_from_rgb_cube(data)

        # NAXIS3(カラー軸)を除去し、2次元画像として書き出す
        if "NAXIS3" in header:
            del header["NAXIS3"]
        header["NAXIS"] = 2
        header["HISTORY"] = "Green plane extracted from pre-debayered RGB cube"

    else:
        print(
            f"エラー: 想定外のデータ次元です: shape={data.shape}",
            file=sys.stderr,
        )
        return 1

    fits.writeto(args.output, green, header, overwrite=True)

    h, w = green.shape[-2], green.shape[-1]
    print(f"書き出し完了: {args.output} ({w}x{h})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
