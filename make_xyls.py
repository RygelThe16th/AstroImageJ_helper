#!/usr/bin/env python3
"""
make_xyls.py

平文の「X Y」座標リスト(1行1星、空白区切り)を、
solve-fieldが読めるFITSバイナリテーブル(xylist)に変換する。

使い方:
    python3 make_xyls.py input.txt output.xyls
"""
import sys

import numpy as np
from astropy.io import fits


def main():
    if len(sys.argv) != 3:
        print("usage: make_xyls.py input.txt output.xyls", file=sys.stderr)
        return 1

    txt_path, out_path = sys.argv[1], sys.argv[2]

    xs, ys = [], []
    with open(txt_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            xs.append(float(parts[0]))
            ys.append(float(parts[1]))

    if not xs:
        print("エラー: 座標が1つも読めませんでした", file=sys.stderr)
        return 1

    col_x = fits.Column(name="X", format="D", array=np.array(xs))
    col_y = fits.Column(name="Y", format="D", array=np.array(ys))
    table_hdu = fits.BinTableHDU.from_columns([col_x, col_y])
    hdul = fits.HDUList([fits.PrimaryHDU(), table_hdu])
    hdul.writeto(out_path, overwrite=True)
    print(f"{len(xs)} 個の星を{out_path}に書き出しました")
    return 0


if __name__ == "__main__":
    sys.exit(main())
