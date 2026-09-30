#!/usr/bin/env python3
"""
fit_minimum.py

tbl_to_csv.rbの出力(datetime_jst, mag, plot_y の3列CSV)から、
食変光星の極小時刻を放物線フィットで求める。

やり方:
    1. 極小候補(生データの中で最も暗い=magが最大の点)を仮の中心とする。
    2. その前後 --window-minutes 分以内の点だけを使って
       mag = a*(t - t0)^2 + b の放物線を最小二乗フィットする
       (tは秒単位の相対時刻)。
    3. フィットの頂点(dmag/dt=0となるt0)を極小時刻とする。
       t0の不確かさは、フィットの共分散行列から求める。

放物線はあくまで極小"付近"での近似(実際の食の光度曲線は、深い場所から
離れるとほぼ直線的な増光・減光になり、放物線からは外れていく)。
--window-minutesを広げすぎると、直線的な部分まで含めてしまい、
かえって極小時刻の推定がずれるので注意。診断用プロット(fit_minimum.png)で、
フィット範囲内の残差が体系的に曲がっていないか(=窓が広すぎないか)を
必ず目で確認すること。

使い方:
    python3 fit_minimum.py rzcas.csv
    python3 fit_minimum.py rzcas.csv --window-minutes 15
    python3 fit_minimum.py rzcas.csv --center "2024-12-06 20:30:00" --window-minutes 20

出力:
    標準出力に、極小時刻(JST)とJD、その1シグマ不確かさを表示する。
    fit_minimum.png に、使ったデータ・フィット曲線・残差を描画する
    (--output-plotで保存先を変更可能)。
"""

import argparse
import sys
from datetime import datetime, timedelta

import numpy as np


def load_csv(path, min_mag=None, max_mag=None):
    times, mags = [], []
    with open(path, encoding='utf-8') as f:
        header = f.readline().strip().split(',')
        dt_idx = header.index('datetime_jst') if 'datetime_jst' in header \
            else header.index('JD')
        mag_idx = header.index('mag')
        is_jd = header[dt_idx] == 'JD'
        n_dropped = 0
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(',')
            mag = float(parts[mag_idx])
            # 基本的な健全性チェック(常時適用): mag<=0やmag>=30は、測定失敗
            # (センタリング失敗、ゼロ除算等)による非物理的な値である可能性が
            # 非常に高いので、フィットにもプロットにも使わない。これを
            # 混入させたまま自動スケールのプロットを作ると、Y軸が壊れて
            # 肝心のデータが見えなくなることを実測で確認したため。
            if not (0 < mag < 30):
                n_dropped += 1
                continue
            if min_mag is not None and mag < min_mag:
                n_dropped += 1
                continue
            if max_mag is not None and mag > max_mag:
                n_dropped += 1
                continue
            if is_jd:
                # JD(ユリウス日)のまま。秒への換算は呼び出し側で行う。
                times.append(float(parts[dt_idx]))
            else:
                times.append(datetime.strptime(parts[dt_idx],
                                                '%Y-%m-%d %H:%M:%S'))
            mags.append(mag)
    if n_dropped:
        print(f'{n_dropped}行を、非物理的な等級値として除外しました '
              f'(0<mag<30の範囲外、または--min-mag/--max-magの範囲外)。')
    return times, np.array(mags), is_jd


def to_seconds(times, ref):
    return np.array([(t - ref).total_seconds() for t in times])


def jd_from_datetime_jst(dt):
    """JST(UTC+9)のdatetimeからユリウス日を計算する(概算、うるう秒等は無視)。"""
    utc = dt - timedelta(hours=9)
    a = (14 - utc.month) // 12
    y = utc.year + 4800 - a
    m = utc.month + 12 * a - 3
    jdn = (utc.day + (153 * m + 2) // 5 + 365 * y + y // 4
           - y // 100 + y // 400 - 32045)
    frac = (utc.hour - 12) / 24 + utc.minute / 1440 + utc.second / 86400
    return jdn + frac


def fit_parabola(t_sec, mag):
    """mag = a*t^2 + b*t + c を最小二乗フィットし、
    頂点位置t0とその1シグマ不確かさを返す。

    注意: 等級は数値が大きいほど暗いので、"光度の極小(一番暗い瞬間)"は
    mag-t平面では山(上に凸、a<0)になる。下に凸(a>0)になった場合は、
    窓の中に本物の極小が無い(例えば単調増光・減光している範囲だけを
    切り取ってしまった)ことを疑うべき。"""
    coeffs, cov = np.polyfit(t_sec, mag, 2, cov=True)
    a, b, _c = coeffs
    if a >= 0:
        raise ValueError(
            '放物線が上に凸(極小=等級のピーク)になりませんでした。窓の中に'
            '極小が含まれていないか、点数が少なすぎる可能性があります。'
            '--window-minutesや--centerを調整してください。')
    t0 = -b / (2 * a)

    # t0 = -b/(2a) の誤差伝播(a, bの共分散を考慮)
    var_a, var_b = cov[0, 0], cov[1, 1]
    cov_ab = cov[0, 1]
    dt0_da = b / (2 * a**2)
    dt0_db = -1 / (2 * a)
    var_t0 = (dt0_da**2 * var_a + dt0_db**2 * var_b
              + 2 * dt0_da * dt0_db * cov_ab)
    sigma_t0 = np.sqrt(var_t0) if var_t0 > 0 else float('nan')
    return t0, sigma_t0, coeffs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('csv_path', help='tbl_to_csv.rbの出力CSV')
    ap.add_argument('--window-minutes', type=float, default=20.0,
                     help='極小候補の前後、フィットに使う時間窓(分, default: 20)')
    ap.add_argument('--center', default=None,
                     help='フィットの中心時刻を手動指定する(JST, '
                          '"YYYY-MM-DD HH:MM:SS"形式)。省略時は生データ中'
                          '最も暗い点を中心にする。')
    ap.add_argument('--output-plot', default='fit_minimum.png',
                     help='診断用プロットの保存先 (default: fit_minimum.png)')
    ap.add_argument('--target-name', default=None,
                     help='グラフのタイトルに使う対象天体名(例: "VW Cep")。'
                          '省略時は"Minimum-time fit"という汎用タイトルになる。')
    ap.add_argument('--min-mag', type=float, default=None,
                     help='これより明るい(数値が小さい)行は除外する'
                          '(0<mag<30の基本チェックに加えて適用)。')
    ap.add_argument('--max-mag', type=float, default=None,
                     help='これより暗い(数値が大きい)行は除外する'
                          '(0<mag<30の基本チェックに加えて適用)。')
    args = ap.parse_args()

    times, mags, is_jd = load_csv(args.csv_path, min_mag=args.min_mag,
                                   max_mag=args.max_mag)
    if is_jd:
        print('エラー: このスクリプトは現状 datetime_jst 列を前提にしています。'
              'tbl_to_csv.rb --jst で出力したCSVを使ってください。',
              file=sys.stderr)
        return 1

    if args.center:
        center_dt = datetime.strptime(args.center, '%Y-%m-%d %H:%M:%S')
    else:
        deepest_idx = int(np.argmax(mags))
        center_dt = times[deepest_idx]
        print(f'中心時刻を自動選択(生データ中で最も暗い点): '
              f'{center_dt} (mag={mags[deepest_idx]:.4f})')

    window_sec = args.window_minutes * 60.0
    t_all_sec = np.array([(t - center_dt).total_seconds() for t in times])
    mask = np.abs(t_all_sec) <= window_sec
    n_used = mask.sum()
    if n_used < 5:
        print(f'エラー: 窓内の点が{n_used}個しかありません(最低5点必要)。'
              '--window-minutesを広げてください。', file=sys.stderr)
        return 1

    t_fit = t_all_sec[mask]
    mag_fit = mags[mask]

    t0_sec, sigma_t0_sec, coeffs = fit_parabola(t_fit, mag_fit)
    t0_dt = center_dt + timedelta(seconds=t0_sec)
    jd_min = jd_from_datetime_jst(t0_dt)
    sigma_jd = sigma_t0_sec / 86400.0

    print()
    print(f'使用した点数: {n_used} / {len(mags)} '
          f'(中心 ±{args.window_minutes:.1f}分の窓)')
    print(f'極小時刻(JST): {t0_dt.strftime("%Y-%m-%d %H:%M:%S")} '
          f'± {sigma_t0_sec:.0f}秒 ({sigma_t0_sec/60:.2f}分)')
    print(f'極小時刻(JD, HJD/BJD補正なしの素のJD): '
          f'{jd_min:.6f} ± {sigma_jd:.6f}')
    print()
    print('注意: これは地球時系のJD(HJD/BJDではない)です。O-C計算等で'
          '他の観測と正確に比較する場合は、光行時間補正(HJD/BJD変換)を'
          '別途行ってください。')

    # 診断用プロット
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(8, 7), sharex=True,
            gridspec_kw={'height_ratios': [3, 1]})

        star_label = args.target_name if args.target_name else None
        title = f'{star_label} minimum-time fit' if star_label \
            else 'Minimum-time fit'

        t_all_min = t_all_sec / 60.0
        ax1.scatter(t_all_min[~mask], mags[~mask], s=12, c='lightgray',
                    label='excluded (outside window)')
        ax1.scatter(t_all_min[mask], mags[mask], s=16, c='#1f4e99',
                    label='used for fit')
        t_curve = np.linspace(t_fit.min(), t_fit.max(), 400)
        mag_curve = np.polyval(coeffs, t_curve)
        ax1.plot(t_curve / 60.0, mag_curve, c='crimson', lw=1.5,
                  label='parabola fit')
        ax1.axvline(t0_sec / 60.0, c='crimson', ls='--', lw=1,
                    label=f'minimum (absolute clock time) '
                          f't0={t0_dt.strftime("%H:%M:%S")}')
        ax1.invert_yaxis()
        ax1.set_ylabel('V magnitude')
        ax1.set_title(title)
        ax1.legend(fontsize=8)
        ax1.grid(alpha=0.3)

        resid = mag_fit - np.polyval(coeffs, t_fit)
        ax2.scatter(t_fit / 60.0, resid, s=14, c='#1f4e99')
        ax2.axhline(0, c='gray', lw=1)
        ax2.set_xlabel(f'minutes from center (center={center_dt.strftime("%H:%M:%S")} '
                       f'clock time; 0 here is NOT the fitted minimum)',
                       fontsize=9)
        ax2.set_ylabel('residual (mag)')
        ax2.grid(alpha=0.3)

        fig.tight_layout()
        fig.savefig(args.output_plot, dpi=130)
        print(f'診断用プロットを書き出しました: {args.output_plot}')
        print('残差が中心付近で系統的に曲がっていないか確認してください'
              '(曲がっていれば--window-minutesを狭める、'
              '端で暴れていれば広げるか除外を検討)。')
    except ImportError:
        print('matplotlibが無いためプロットは省略します。', file=sys.stderr)

    return 0


if __name__ == '__main__':
    sys.exit(main())
