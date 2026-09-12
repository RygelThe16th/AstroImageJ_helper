#!/usr/bin/env python3
"""
stack_and_solve.py

441枚のような大量の個別ライトフレーム(1枚あたりでは星が少なすぎて
plate solveできないことがある)を、連続するNフレームずつシグマクリップ
平均でスタックしてSN比を稼ぎ、そのスタック画像を「本物のFITS画像」として
solve-fieldに直接渡してWCSを解く。

これによりAIJの星検出(Noise Tol等)を一切経由せず、
また明るさ順のソートの問題も気にする必要がなくなる
(solve-field自身が画像から検出するため)。

解けたWCSは、そのグループに属する元の個別フレーム全部の
ヘッダーにコピーして、別ディレクトリに書き出す
(元ファイルは一切変更しない)。

あるグループが解けなかった場合は、隣接する次のグループと自動的に
合体して(枚数を増やして)再挑戦する。これを--max-group-sizeで
指定した上限枚数まで繰り返し、それでも解けなければ最終的に
失敗として扱う(WCSを付けずにコピーのみ行う)。

前提:
    - astrometry.net (solve-field) がインストール済みで、
      該当スケールのインデックスファイルも配置済みであること
    - pip install astropy numpy --break-system-packages

使い方:
    python3 stack_and_solve.py \
        --input-dir /path/to/lights \
        --output-dir /path/to/lights_with_wcs \
        --pattern "Light_rz_cas_*.fit*" \
        --group-size 5 \
        --scale-low 3.74 --scale-high 5.74 --scale-units arcsecperpix \
        --backend-config /etc/astrometry.cfg

出力:
    --output-dir 以下に、入力と同じファイル名で、WCSヘッダーが
    追加(上書き)された個別フレームが書き出される。
    最終的に解決できなかったグループのフレームは、WCSを付けずに
    そのままコピーされる(一覧は最後にまとめて表示される)。
    スタック中間ファイルやsolve-fieldの作業ファイルは
    --work-dir (デフォルトは一時ディレクトリ) に残る。
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import deque

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip
from scipy.ndimage import shift as ndi_shift
from skimage.registration import phase_cross_correlation

# solve-fieldの出力等でヘッダーに残っている可能性のある、WCS関連キーワードの
# プレフィックス/完全一致名。スタック画像を作る前、および結果を書き戻す前の
# 両方で、これらを一旦除去してから新しいWCSを書き込む。
WCS_KEY_PREFIXES = ('CRVAL', 'CRPIX', 'CD1_', 'CD2_', 'CDELT', 'CTYPE',
                    'CUNIT', 'PV1_', 'PV2_', 'PC1_', 'PC2_',
                    'A_', 'B_', 'AP_', 'BP_')
WCS_KEY_EXACT = ('WCSAXES', 'RADESYS', 'EQUINOX', 'LONPOLE', 'LATPOLE')


def strip_wcs_keywords(header):
    for key in list(header.keys()):
        if key.startswith(WCS_KEY_PREFIXES) or key in WCS_KEY_EXACT:
            del header[key]


def find_frames(input_dir, pattern):
    paths = sorted(glob.glob(os.path.join(input_dir, pattern)))
    if not paths:
        print(f"エラー: {input_dir} 内に '{pattern}' に一致するファイルがありません。",
              file=sys.stderr)
        sys.exit(1)
    return paths


def group_frames(paths, group_size):
    """連続するgroup_size枚ずつのグループに分割する。最後の端数グループも
    そのまま(枚数が少なくても)処理する。"""
    return [paths[i:i + group_size] for i in range(0, len(paths), group_size)]


def load_group_data(paths):
    """グループ内の各フレームの主HDUデータを読み込み、(N,H,W)のarrayで返す。
    最初のフレームのヘッダーも代表として返す(スタック画像のNAXIS等に使う)。"""
    stack = []
    header = None
    for p in paths:
        with fits.open(p, memmap=False) as hdul:
            data = hdul[0].data
            if header is None:
                header = hdul[0].header.copy()
            stack.append(data.astype(np.float64))
    return np.stack(stack, axis=0), header


def pick_reference_index(data_cube):
    """グループ内で最も"良い"(コントラストが高い=星がはっきり写っている
    可能性が高い)フレームを基準フレームとして選ぶ。単純に先頭フレームを
    固定で使うと、そのフレームがたまたま条件の悪い(星がほとんど見えない)
    ものだった場合に、残り全部の位置合わせが道連れで破綻することが
    実測で疑われたため。"背景を除いた最大値"が最も高いフレームを選ぶ。"""
    scores = []
    for frame in data_cube:
        bg = np.nanmedian(frame)
        scores.append(np.nanmax(frame) - bg)
    return int(np.argmax(scores))


def correlation_quality(reference, shifted_image, shift):
    """shift適用後の画像が、referenceとどれだけ実際に重なっているかを
    正規化相互相関係数(-1〜1、1に近いほど良く重なる)で直接測る。

    「シフト量が時間に対してなめらかか」という間接的な前提(以前試した
    方式)は、望遠鏡が意図的に行う"ディザリング"(数フレームおきの、
    小さくない意図的な視野移動)がある場合に、本物の大きな移動を
    誤検出として弾いてしまうことが分かった。この指標は移動量の大小に
    関わらず、シフト後に実際に画像同士が一致するかどうかだけを見るので、
    本物の大きなジャンプと、単なる相関の誤検出(ノイズ同士がたまたま
    ある位置で"一致した"ことにされてしまう)を正しく区別できる。"""
    dy, dx = shift
    pad_y = int(np.ceil(abs(dy))) + 2
    pad_x = int(np.ceil(abs(dx))) + 2
    h, w = reference.shape
    ys, ye, xs, xe = pad_y, h - pad_y, pad_x, w - pad_x
    if ye <= ys or xe <= xs:
        return 0.0
    a = reference[ys:ye, xs:xe] - np.mean(reference[ys:ye, xs:xe])
    b = shifted_image[ys:ye, xs:xe] - np.mean(shifted_image[ys:ye, xs:xe])
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    if denom == 0:
        return 0.0
    return float(np.sum(a * b) / denom)


def fill_unreliable_shifts_from_neighbors(shifts_used, reliable):
    """信頼できると判定された(quality検査を通った)フレームの値だけを
    使って、信頼できないフレームのシフト値を前後から線形補間で埋める。
    信頼できるフレームが1つも無い場合は補間できないので、そのまま返す。"""
    n = len(shifts_used)
    reliable_idx = [i for i in range(n) if reliable[i]]
    if not reliable_idx:
        return shifts_used
    idx = np.arange(n)
    dy = np.array([shifts_used[i][0] for i in range(n)])
    dx = np.array([shifts_used[i][1] for i in range(n)])
    dy_filled = np.interp(idx, [i for i in reliable_idx],
                           [dy[i] for i in reliable_idx])
    dx_filled = np.interp(idx, [i for i in reliable_idx],
                           [dx[i] for i in reliable_idx])
    return [(float(dy_filled[i]) if not reliable[i] else shifts_used[i][0],
             float(dx_filled[i]) if not reliable[i] else shifts_used[i][1])
            for i in range(n)]


def align_stack(data_cube, upsample_factor=10, max_shift=None,
                 quality_threshold=0.1, verbose=False):
    """(N,H,W)のデータキューブの各フレームを、グループ内で最も良い
    (コントラストが高い)フレームを基準として、サブピクセル位置合わせする
    (位相相関法)。

    フレーム間で星の位置がずれたまま単純に平均・スタックすると、星の光が
    複数ピクセルに"にじんで薄まり"、スタック枚数を増やすほど逆に検出しにくく
    なることが実測で確認されている(このデータはレジストレーション未実施の
    生フレームであるため)。位置合わせしてからスタックすることで、
    枚数を増やすほど本来のSN比向上が得られるようにする。

    基準フレームから離れたフレームを直接1回の相関で合わせようとすると、
    星が少なく暗いフレームほど不安定になることを実測で確認したため、
    「これまでに位置合わせ済みのフレームを足し合わせた、育っていく
    合成画像」を相手に、外側へ向かって順番に相関を取る方式にしている。

    各ステップの信頼度は、シフト量の大小や"なめらかさ"ではなく、
    correlation_quality()で実際にどれだけ画像同士が重なるかを直接測って
    判定する(quality_threshold未満は信頼できないとみなす)。信頼できない
    フレームは、育っていく合成画像の更新には使わない(誤った位置合わせで
    以降のフレームの基準そのものを汚染しないため)が、そのフレーム自身の
    最終的なシフト値は、信頼できる前後のフレームから補間して埋める。"""
    ref_idx = pick_reference_index(data_cube)
    n = data_cube.shape[0]
    shifts_used = [None] * n
    quality_used = [None] * n
    shifts_used[ref_idx] = (0.0, 0.0)
    quality_used[ref_idx] = 1.0

    def estimate_shift(template, image):
        try:
            shift_yx, _error, _diffphase = phase_cross_correlation(
                template, image, upsample_factor=upsample_factor)
        except Exception:
            shift_yx = np.array([0.0, 0.0])
        if max_shift is not None and (abs(shift_yx[0]) > max_shift
                                       or abs(shift_yx[1]) > max_shift):
            shift_yx = np.array([0.0, 0.0])
        shifted = ndi_shift(image, shift_yx, order=3, mode='nearest')
        quality = correlation_quality(template, shifted, shift_yx)
        return shift_yx, quality

    def grow_direction(indices):
        running_stack = data_cube[ref_idx].astype(np.float64).copy()
        running_count = 1
        for i in indices:
            template = running_stack / running_count
            step, quality = estimate_shift(template, data_cube[i])
            shifts_used[i] = (float(step[0]), float(step[1]))
            quality_used[i] = quality
            if quality >= quality_threshold:
                aligned_i = ndi_shift(data_cube[i], step, order=3,
                                       mode='nearest')
                running_stack += aligned_i
                running_count += 1
            # qualityが低い場合は、育っていく合成画像を汚染しないよう
            # このフレームは足し込まず、次のフレームも同じtemplateで試す。

    grow_direction(range(ref_idx + 1, n))
    grow_direction(range(ref_idx - 1, -1, -1))

    reliable = [q >= quality_threshold for q in quality_used]
    reliable[ref_idx] = True
    shifts_used = fill_unreliable_shifts_from_neighbors(shifts_used, reliable)

    aligned = np.empty_like(data_cube)
    for i in range(n):
        dy, dx = shifts_used[i]
        if dy == 0.0 and dx == 0.0:
            aligned[i] = data_cube[i]
        else:
            aligned[i] = ndi_shift(data_cube[i], (dy, dx), order=3,
                                    mode='nearest')
    if verbose:
        print(f"    基準フレーム: {ref_idx}番目")
        for i, s in enumerate(shifts_used):
            q = quality_used[i]
            flag = '' if reliable[i] else f' (信頼度低, quality={q:.3f})'
            print(f"    frame[{i}] shift=({s[0]:+.2f}, {s[1]:+.2f}){flag}")
    return aligned, shifts_used, reliable


def sigma_clipped_stack(data_cube, sigma=3.0):
    """(N,H,W)のデータキューブを、各ピクセルごとにNフレーム方向へ
    シグマクリップした上で平均する。

    cenfunc/stdfuncは既定(平均/標準偏差)ではなく、外れ値自体に引きずられない
    median/mad_std(中央絶対偏差ベース)を使う。既定のままだと、宇宙線などの
    単一の極端な外れ値がstd自体を吊り上げてしまい、標本数が少ない(数枚)場合に
    全く弾かれないことがある(実測で確認済み)。"""
    clipped = sigma_clip(data_cube, sigma=sigma, axis=0, masked=True,
                          cenfunc='median', stdfunc='mad_std')
    combined = np.ma.mean(clipped, axis=0)
    return np.ma.filled(combined, np.nan)


def write_stack_fits(path, data, header):
    out_header = header.copy()
    strip_wcs_keywords(out_header)
    fits.writeto(path, data.astype(np.float32), out_header, overwrite=True)


def run_solve_field(image_path, work_dir, scale_low, scale_high, scale_units,
                     backend_config, tweak_order, extra_args):
    args = [
        'solve-field', '--overwrite', '--no-plots',
        '-D', work_dir,
        '--backend-config', backend_config,
        '-u', scale_units, '-L', str(scale_low), '-H', str(scale_high),
        '--tweak-order', str(tweak_order),
    ]
    args += extra_args
    args.append(image_path)
    started = time.time()
    proc = subprocess.run(args, capture_output=True, text=True)
    elapsed = time.time() - started
    return proc, elapsed


def load_solved_wcs_header(work_dir, stem):
    wcs_path = os.path.join(work_dir, f'{stem}.wcs')
    if not os.path.exists(wcs_path):
        return None
    with fits.open(wcs_path, memmap=False) as hdul:
        return hdul[0].header.copy()


def apply_wcs_to_frame(src_path, dst_path, wcs_header, frame_shift=None):
    """wcs_headerをsrc_pathのフレームに適用してdst_pathへ書き出す。

    frame_shift(dy, dx)は、このフレームがスタック前のアライメントで
    基準フレームに対して"ずらされた"量。wcs_headerはズラした後の
    (=基準フレームの)座標系で解かれているので、それをズラす前の
    このフレーム自身の生ピクセル座標系に正しく対応させるには、
    CRPIX1/CRPIX2をこのシフト量ぶんだけ補正する必要がある
    (これを怠ると、シフト量が大きいフレームほどWCSが実際の星の位置から
    ズレ、Multi-Apertureの"No signal for centroid"エラーの原因になる
    ことを確認した)。frame_shiftがNone、または(0,0)なら補正不要。"""
    with fits.open(src_path, memmap=False) as hdul:
        data = hdul[0].data
        header = hdul[0].header.copy()
    strip_wcs_keywords(header)
    for card in wcs_header.cards:
        if card.keyword in ('SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2',
                             'EXTEND', 'COMMENT', 'HISTORY', ''):
            continue
        header[card.keyword] = (card.value, card.comment)
    if frame_shift is not None:
        dy, dx = frame_shift
        if dy != 0.0 and 'CRPIX2' in header:
            header['CRPIX2'] = header['CRPIX2'] - dy
        if dx != 0.0 and 'CRPIX1' in header:
            header['CRPIX1'] = header['CRPIX1'] - dx
    fits.writeto(dst_path, data, header, overwrite=True)


def solve_group(group, group_dir, args):
    """1グループ(フレームのパスのリスト)を(必要なら位置合わせしてから)
    スタックしてsolve-fieldにかける。
    (解けたWCSヘッダー(astropy Header)またはNone, 所要秒数,
     groupと同じ順番のシフト量[(dy,dx),...]またはNone,
     同じ順番の信頼度[bool,...]またはNone)を返す。
     信頼度がFalseのフレームは、シフト量自体は(なめらかさを仮定した
     推定値で)埋まっているが、実際の相関では大きく外れていたフレームで
     あることを示す。呼び出し側は、そのようなフレームには個別のWCSを
     適用しない(付けずにコピーのみ行う)べき。"""
    os.makedirs(group_dir, exist_ok=True)
    stack_path = os.path.join(group_dir, 'stack.fits')

    data_cube, header = load_group_data(group)
    shifts = None
    reliable = None
    if not args.no_align and data_cube.shape[0] > 1:
        data_cube, shifts, reliable = align_stack(
            data_cube, upsample_factor=args.align_upsample_factor,
            max_shift=args.align_max_shift,
            quality_threshold=args.align_quality_threshold,
            verbose=args.verbose_align)
        max_abs_shift = max(max(abs(dy), abs(dx)) for dy, dx in shifts)
        if max_abs_shift > 0.5:
            print(f"  位置合わせ: 最大シフト量 {max_abs_shift:.2f}px "
                  f"({len(group)}枚)")
    stacked = sigma_clipped_stack(data_cube, sigma=args.sigma)
    write_stack_fits(stack_path, stacked, header)

    _, elapsed = run_solve_field(
        stack_path, group_dir, args.scale_low, args.scale_high,
        args.scale_units, args.backend_config, args.tweak_order,
        args.extra_solve_field_arg)

    wcs_header = load_solved_wcs_header(group_dir, 'stack')
    return wcs_header, elapsed, shifts, reliable


def main():
    # 標準出力がtee等のパイプに繋がれると、Pythonはデフォルトでフルバッファリング
    # (端末に直結時の行バッファリングと違い、数KB溜まるかプログラム終了まで
    # 出力が出てこない)になる。ここで明示的に行バッファリングへ切り替えて、
    # tee経由でもリアルタイムに1行ずつログが見えるようにする。
    sys.stdout.reconfigure(line_buffering=True)

    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--input-dir', required=True)
    ap.add_argument('--output-dir', required=True)
    ap.add_argument('--pattern', default='*.fits')
    ap.add_argument('--group-size', type=int, default=5)
    ap.add_argument('--scale-low', type=float, default=3.74)
    ap.add_argument('--scale-high', type=float, default=5.74)
    ap.add_argument('--scale-units', default='arcsecperpix')
    ap.add_argument('--backend-config', default='/etc/astrometry.cfg')
    ap.add_argument('--tweak-order', type=int, default=2)
    ap.add_argument('--sigma', type=float, default=3.0,
                     help='スタック時のシグマクリップ閾値 (default: 3.0)')
    ap.add_argument('--work-dir', default=None,
                     help='solve-fieldの作業ディレクトリ(未指定なら一時ディレクトリ)')
    ap.add_argument('--extra-solve-field-arg', action='append', default=[],
                     help='solve-fieldにそのまま追加で渡す引数(複数指定可)')
    ap.add_argument('--max-group-size', type=int, default=None,
                     help='失敗グループを隣接グループと自動合体させて再挑戦する際の'
                          '上限枚数(未指定なら--group-sizeの4倍)。この枚数を'
                          '超えてもなお解けなければ、最終的に失敗として扱う。')
    ap.add_argument('--no-align', action='store_true',
                     help='スタック前の位置合わせ(位相相関法)を無効化する。'
                          '既にレジストレーション済みのデータなど、フレーム間で'
                          '星が動かないことが分かっている場合のみ指定すること。'
                          '通常は有効(デフォルト)のままにする。')
    ap.add_argument('--align-upsample-factor', type=int, default=10,
                     help='位置合わせのサブピクセル精度(大きいほど精密だが遅い。'
                          'default: 10 = 1/10px精度)')
    ap.add_argument('--align-max-shift', type=float, default=None,
                     help='位置合わせで許容する最大シフト量(px)の安全上限。'
                          '指定した場合、これを超える推定値は問答無用でshift=0'
                          'として扱う。通常は指定不要(未指定=上限無し)。'
                          '信頼度判定は--align-quality-thresholdの方で行われる'
                          '(shift量の大小そのものは信頼度の判断材料にしない。'
                          '望遠鏡が意図的に行う大きなディザリングを、誤検出と'
                          '誤って弾いてしまうことがあったため)。')
    ap.add_argument('--align-quality-threshold', type=float, default=0.1,
                     help='位置合わせの信頼度判定に使う、正規化相互相関係数の'
                          '下限(default: 0.1)。シフト適用後に実際にどれだけ'
                          '画像同士が重なるかを直接測る指標で、これを下回ると'
                          'そのフレームの位置合わせは信頼できないと判定される。')
    ap.add_argument('--skip', type=int, default=0,
                     help='ファイル名順に並べた先頭から、この枚数だけスキップして'
                          'から処理を始める(default: 0)。特定区間だけを'
                          '切り出してテストしたい場合に使う。')
    ap.add_argument('--limit', type=int, default=None,
                     help='--skip適用後、最大この枚数までしか処理しない'
                          '(未指定なら残り全部)。')
    ap.add_argument('--verbose-align', action='store_true',
                     help='各グループの位置合わせで、基準フレームと'
                          '各フレームごとのシフト量を全部表示する'
                          '(デバッグ用。通常は最大値だけ表示される)。')
    ap.add_argument('--write-unsolved', action='store_true',
                     help='最終的に解決できなかったグループのフレームも、'
                          'WCS無しでそのまま出力ディレクトリにコピーする。'
                          'デフォルトでは書き出さない(出力ディレクトリには'
                          '解決できたフレームだけが残る)。')
    args = ap.parse_args()
    max_group_size = args.max_group_size or args.group_size * 4

    os.makedirs(args.output_dir, exist_ok=True)
    work_dir_root = args.work_dir or tempfile.mkdtemp(prefix='stack_and_solve_')
    os.makedirs(work_dir_root, exist_ok=True)

    frame_paths = find_frames(args.input_dir, args.pattern)
    if args.skip or args.limit is not None:
        end = None if args.limit is None else args.skip + args.limit
        before = len(frame_paths)
        frame_paths = frame_paths[args.skip:end]
        print(f"--skip {args.skip} --limit {args.limit} を適用: "
              f"{before}枚中 {len(frame_paths)}枚を対象にする "
              f"({frame_paths[0] if frame_paths else '(該当なし)'} 〜 "
              f"{frame_paths[-1] if frame_paths else ''})")
        if not frame_paths:
            print("エラー: --skip/--limit の結果、対象フレームが0枚になりました。",
                  file=sys.stderr)
            sys.exit(1)
    initial_groups = group_frames(frame_paths, args.group_size)
    print(f"{len(frame_paths)}枚を{len(initial_groups)}グループ"
          f"(1グループ{args.group_size}枚, 失敗時は最大{max_group_size}枚まで"
          f"隣と合体して再挑戦)に分割")
    print(f"作業ディレクトリ: {work_dir_root}")

    # 隣接グループと合体しながら処理するため、キュー(deque)で扱う。
    # 各要素は(表示・作業ディレクトリ名に使うラベル, フレームパスのリスト)。
    queue = deque((f'{i:04d}', g) for i, g in enumerate(initial_groups))

    n_ok, n_fail = 0, 0
    n_unreliable_frames = 0
    unreliable_frame_names = []
    failed_groups = []
    attempt_seq = 0

    while queue:
        label, group = queue.popleft()
        attempt_seq += 1
        group_dir = os.path.join(work_dir_root, f'group_{label}')
        wcs_header, elapsed, shifts, reliable = solve_group(group, group_dir, args)
        names = [os.path.basename(p) for p in group]

        if wcs_header is not None:
            n_ok += 1
            print(f"[group {label}] 解決成功 ({elapsed:.2f}秒, {len(group)}枚: "
                  f"{names[0]} 〜 {names[-1]})")
            for idx, p in enumerate(group):
                dst = os.path.join(args.output_dir, os.path.basename(p))
                is_reliable = reliable[idx] if reliable is not None else True
                if not is_reliable:
                    # 位置合わせの信頼度は低いが、AIJのMulti-Apertureは
                    # WCS無しのフレームに出会うとスキップではなく処理全体を
                    # 中断してしまう(実機で確認済み)ため、WCSを付けないより
                    # 多少不正確でも(なめらかなドリフトを仮定した推定値で)
                    # 付けておく方が実用上マシと判断し、必ずWCSを書き込む。
                    # どのフレームが低信頼度だったかは記録だけ残す。
                    n_unreliable_frames += 1
                    unreliable_frame_names.append(os.path.basename(p))
                frame_shift = shifts[idx] if shifts is not None else None
                apply_wcs_to_frame(p, dst, wcs_header, frame_shift=frame_shift)
            continue

        # 失敗:次のグループが控えていて、合体後も上限枚数以内なら合体して再挑戦
        if queue and len(group) + len(queue[0][1]) <= max_group_size:
            next_label, next_group = queue.popleft()
            merged_label = f'{label}+{next_label}'
            merged_group = group + next_group
            print(f"[group {label}] 解決失敗 ({elapsed:.2f}秒, {len(group)}枚)"
                  f" → 次のグループ({next_label})と合体して再挑戦"
                  f"({len(merged_group)}枚)")
            queue.appendleft((merged_label, merged_group))
            continue

        # これ以上合体できない(上限到達、またはもう後続グループが無い)ので
        # 最終的な失敗として確定する
        n_fail += 1
        failed_groups.append((label, names))
        print(f"[group {label}] 解決失敗 ({elapsed:.2f}秒, {len(group)}枚: "
              f"{names[0]} 〜 {names[-1]}) → これ以上合体できないため断念",
              file=sys.stderr)
        if args.write_unsolved:
            for p in group:
                dst = os.path.join(args.output_dir, os.path.basename(p))
                shutil.copy2(p, dst)

    print()
    print(f"完了: {n_ok}グループ成功 / {n_fail}グループ最終失敗 "
          f"(計{len(frame_paths)}枚, 試行回数{attempt_seq})")
    if n_unreliable_frames:
        print(f"うち、位置合わせの信頼度が低いと判定されたが、"
              f"(なめらかなドリフトを仮定した推定値で)WCSは付けたフレーム: "
              f"{n_unreliable_frames}枚")
        for name in unreliable_frame_names:
            print(f"  {name}")
        print("→ これらのフレームは他より測光精度が落ちている可能性があります。"
              "光度曲線で明らかに外れた点が出たら、まずこの一覧を疑ってください。")
    if failed_groups:
        skip_note = "(出力ディレクトリには書き出していません)" if not args.write_unsolved \
            else "(WCS無しでコピー済み)"
        print(f"最終的に解けなかったグループ {skip_note}:")
        for label, names in failed_groups:
            print(f"  group {label}: {names[0]} 〜 {names[-1]} ({len(names)}枚)")
        print(f"→ --max-group-size {max_group_size} でも解けなかった区間です。"
              "露光時間や視野の星密度そのものの限界の可能性があるので、"
              "該当区間だけさらに広い範囲で手動合体するか、個別に確認してください。"
              "(--write-unsolvedを付けるとWCS無しでもコピーされます)")


if __name__ == '__main__':
    sys.exit(main())
