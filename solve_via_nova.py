#!/usr/bin/env python3
"""
solve_via_nova.py

stack_and_solve.pyが--unreliable-list-outで書き出した「位置合わせの
信頼度が低い」フレーム一覧を読み込み、その各フレームを本物の
nova.astrometry.net(オンラインのastrometry.netサービス)へ個別に
アップロードしてplate solveし、得られたWCSを元のフレームに適用する。

なぜこの少数フレームだけオンライン版を使うのか:
    stack_and_solve.py(ローカルのsolve-field)は、複数フレームを
    スタックしてSN比を稼ぐことが前提の設計になっている。1枚だけの
    フレームをローカルで解こうとすると星の数が足りず失敗しやすい。
    一方、nova.astrometry.netは索引しているインデックスファイルが
    ローカル環境より充実していることが多く、1枚だけでも解けることが
    ある。ただし1枚あたり数秒〜数十秒かかる(混雑状況次第)ため、
    大量のフレームを流す用途には不向き。「低信頼度」フレームは
    通常ごく少数なので、この用途にはちょうど良い。

前提:
    - nova.astrometry.netのアカウントを作り、API Keyを取得しておく
      (https://nova.astrometry.net/api_help を参照)
    - pip install requests astropy

使い方:
    python3 solve_via_nova.py \
        --unreliable-list unreliable.txt \
        --input-dir subsG \
        --output-dir subsG_wcs_fixed \
        --api-key YOUR_API_KEY \
        --scale-low 6.6 --scale-high 8.1 --scale-units arcsecperpix

出力:
    --output-dir 以下に、解けたフレームだけWCSを適用して書き出す。
    解けなかったフレームはコピーもされない(ログに失敗として表示される)。
    元のフレーム(--input-dir側)は一切変更しない。
"""

import argparse
import json
import sys
import time

import requests
from astropy.io import fits

API_BASE = 'https://nova.astrometry.net/api'
# wcs_fileエンドポイントだけ、他の/api/...とは違い/apiプレフィックスが
# 付かない(nova.astrometry.net実サービスの実際の仕様。ruby_ansvr.rbを
# 書いた際にも同じ形で確認済み: curl .../wcs_file/1 に /api は付かない)。
SITE_BASE = 'https://nova.astrometry.net'
WCS_KEY_PREFIXES = ('CRVAL', 'CRPIX', 'CD1_', 'CD2_', 'CDELT', 'CTYPE',
                    'CUNIT', 'PV1_', 'PV2_', 'PC1_', 'PC2_',
                    'A_', 'B_', 'AP_', 'BP_')
WCS_KEY_EXACT = ('WCSAXES', 'RADESYS', 'EQUINOX', 'LONPOLE', 'LATPOLE')


def strip_wcs_keywords(header):
    for key in list(header.keys()):
        if key.startswith(WCS_KEY_PREFIXES) or key in WCS_KEY_EXACT:
            del header[key]


def load_unreliable_list(path):
    """stack_and_solve.pyの--unreliable-list-outの出力形式を読む。
    '#'始まりの行、空行は無視し、残りをファイル名として返す。"""
    names = []
    with open(path, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            names.append(line)
    return names


def request_with_retry(method, url, max_retries=5, base_delay=3.0, **kwargs):
    """公開サービス(nova.astrometry.net)は一時的な503/接続エラーを
    普通に返すことがある。即座に諦めず、指数バックオフで数回まで
    再挑戦する。全部失敗したら最後の例外を送出する。"""
    last_exc = None
    for attempt in range(1, max_retries + 1):
        try:
            resp = method(url, **kwargs)
            resp.raise_for_status()
            return resp
        except requests.exceptions.RequestException as e:
            last_exc = e
            if attempt == max_retries:
                break
            delay = base_delay * (2 ** (attempt - 1))
            print(f'  (一時的なエラー、{delay:.0f}秒後に再挑戦 '
                  f'{attempt}/{max_retries}: {e})', file=sys.stderr)
            time.sleep(delay)
    raise last_exc


def login(api_key):
    resp = request_with_retry(
        requests.post, f'{API_BASE}/login',
        data={'request-json': json.dumps({'apikey': api_key})})
    data = resp.json()
    if data.get('status') != 'success':
        raise RuntimeError(f'ログイン失敗: {data}')
    return data['session']


def upload_file(session, file_path, scale_low, scale_high, scale_units,
                tweak_order):
    opts = {
        'session': session,
        'allow_commercial_use': 'd',
        'allow_modifications': 'd',
        'publicly_visible': 'n',
        'scale_units': scale_units,
        'scale_type': 'ul',
        'scale_lower': scale_low,
        'scale_upper': scale_high,
        'tweak_order': tweak_order,
    }
    with open(file_path, 'rb') as f:
        resp = request_with_retry(
            requests.post, f'{API_BASE}/upload',
            data={'request-json': json.dumps(opts)},
            files={'file': f})
    data = resp.json()
    if data.get('status') != 'success':
        raise RuntimeError(f'アップロード失敗: {data}')
    return data['subid']


def wait_for_job_id(subid, poll_interval, timeout):
    """submissionにjobが割り当てられ、処理が終わるまで待つ。
    (job_id, finished)を返す。job_idがまだ無ければ(None, False)。"""
    started = time.time()
    while time.time() - started < timeout:
        resp = request_with_retry(requests.get,
                                   f'{API_BASE}/submissions/{subid}')
        data = resp.json()
        jobs = data.get('jobs') or []
        if jobs and jobs[0] is not None:
            return jobs[0], bool(data.get('processing_finished'))
        time.sleep(poll_interval)
    return None, False


def wait_for_job_status(job_id, poll_interval, timeout):
    started = time.time()
    while time.time() - started < timeout:
        resp = request_with_retry(requests.get, f'{API_BASE}/jobs/{job_id}')
        status = resp.json().get('status')
        if status in ('success', 'failure'):
            return status
        time.sleep(poll_interval)
    return 'timeout'


def download_wcs_header(job_id):
    resp = request_with_retry(requests.get, f'{SITE_BASE}/wcs_file/{job_id}')
    tmp_path = f'/tmp/solve_via_nova_{job_id}.wcs'
    with open(tmp_path, 'wb') as f:
        f.write(resp.content)
    with fits.open(tmp_path, memmap=False) as hdul:
        return hdul[0].header.copy()


def apply_wcs_to_frame(src_path, dst_path, wcs_header):
    with fits.open(src_path, memmap=False) as hdul:
        data = hdul[0].data
        header = hdul[0].header.copy()
    strip_wcs_keywords(header)
    for card in wcs_header.cards:
        if card.keyword in ('SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2',
                             'EXTEND', 'COMMENT', 'HISTORY', ''):
            continue
        header[card.keyword] = (card.value, card.comment)
    fits.writeto(dst_path, data, header, overwrite=True)


def main():
    import os

    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--unreliable-list', required=True,
                     help='stack_and_solve.pyの--unreliable-list-outの出力ファイル')
    ap.add_argument('--input-dir', required=True,
                     help='元のフレームがあるディレクトリ')
    ap.add_argument('--output-dir', required=True,
                     help='解けたフレームの書き出し先')
    ap.add_argument('--api-key', required=True,
                     help='nova.astrometry.netのAPI Key')
    ap.add_argument('--scale-low', type=float, required=True)
    ap.add_argument('--scale-high', type=float, required=True)
    ap.add_argument('--scale-units', default='arcsecperpix')
    ap.add_argument('--tweak-order', type=int, default=2)
    ap.add_argument('--poll-interval', type=float, default=5.0,
                     help='submission/jobの状態確認の間隔(秒, default: 5.0)')
    ap.add_argument('--timeout', type=float, default=300.0,
                     help='1フレームあたりの最大待ち時間(秒, default: 300.0)')
    args = ap.parse_args()

    sys.stdout.reconfigure(line_buffering=True)

    names = load_unreliable_list(args.unreliable_list)
    if not names:
        print('エラー: 一覧が空です。', file=sys.stderr)
        return 1
    print(f'{len(names)}枚を対象にします: {names}')

    os.makedirs(args.output_dir, exist_ok=True)

    print('nova.astrometry.netにログイン中...')
    session = login(args.api_key)
    print('ログイン成功')

    n_ok, n_fail = 0, 0
    for name in names:
        src_path = os.path.join(args.input_dir, name)
        if not os.path.exists(src_path):
            print(f'[{name}] エラー: 元のファイルが見つかりません: {src_path}',
                  file=sys.stderr)
            n_fail += 1
            continue

        print(f'[{name}] アップロード中...')
        try:
            subid = upload_file(session, src_path, args.scale_low,
                                 args.scale_high, args.scale_units,
                                 args.tweak_order)
        except Exception as e:
            print(f'[{name}] アップロード失敗: {e}', file=sys.stderr)
            n_fail += 1
            continue

        print(f'[{name}] subid={subid}, job割り当て待ち...')
        job_id, _finished = wait_for_job_id(subid, args.poll_interval,
                                             args.timeout)
        if job_id is None:
            print(f'[{name}] タイムアウト: jobが割り当てられませんでした',
                  file=sys.stderr)
            n_fail += 1
            continue

        print(f'[{name}] job_id={job_id}, 解決待ち...')
        status = wait_for_job_status(job_id, args.poll_interval, args.timeout)

        if status != 'success':
            print(f'[{name}] 解決失敗 (status={status})', file=sys.stderr)
            n_fail += 1
            continue

        wcs_header = download_wcs_header(job_id)
        dst_path = os.path.join(args.output_dir, name)
        apply_wcs_to_frame(src_path, dst_path, wcs_header)
        print(f'[{name}] 解決成功、書き出し完了: {dst_path}')
        n_ok += 1

    print()
    print(f'完了: {n_ok}枚成功 / {n_fail}枚失敗 (計{len(names)}枚)')
    return 0 if n_fail == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
