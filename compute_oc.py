#!/usr/bin/env python3
"""
compute_oc.py

観測で求めた食変光星の極小時刻(JST)を、天体の座標を使って
HJD(太陽系重心での光行時間補正込みのユリウス日)に変換し、
既知の元期(Epoch)・周期(Period)から予報時刻との差(O-C)を計算する。

なぜHJD変換が必要か:
    fit_minimum.pyが出す極小時刻は、地球(地心)を基準にした素のJD。
    一方、VSX等に載っている食変光星の元期は通常HJD(太陽系重心が
    基準)。地球の公転により、この2つは天体の方向次第で最大約8分
    (地球軌道の光行時間)もズレるので、正しくO-Cを計算するには
    自分の観測値もHJDに変換してから比較する必要がある。

使い方:
    python3 compute_oc.py \
        --min-time-jst "2024-12-06 20:31:12" \
        --ra "02:48:55.51" --dec "+69:38:03.4" \
        --epoch 2460646.2029 --period 1.19525556

    # fit_minimum.pyの1シグマ不確かさ(秒)も渡すと、O-Cの誤差も計算する
    python3 compute_oc.py \
        --min-time-jst "2024-12-06 20:31:12" --min-time-sigma-sec 17 \
        --ra "02:48:55.51" --dec "+69:38:03.4" \
        --epoch 2460646.2029 --period 1.19525556
"""

import argparse
import sys
from datetime import datetime, timedelta

from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
import astropy.units as u


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--min-time-jst', required=True,
                     help='観測した極小時刻(JST, "YYYY-MM-DD HH:MM:SS"形式)')
    ap.add_argument('--min-time-sigma-sec', type=float, default=None,
                     help='極小時刻の1シグマ不確かさ(秒)。fit_minimum.pyの'
                          '出力にある値をそのまま渡せる。')
    ap.add_argument('--ra', required=True,
                     help='対象天体の赤経(J2000, "HH:MM:SS"形式)')
    ap.add_argument('--dec', required=True,
                     help='対象天体の赤緯(J2000, "+DD:MM:SS"形式)')
    ap.add_argument('--epoch', type=float, required=True,
                     help='元期(HJD)。VSX等に載っている値。')
    ap.add_argument('--period', type=float, required=True,
                     help='周期(日)。')
    ap.add_argument('--tz-offset-hours', type=float, default=9.0,
                     help='入力時刻のタイムゾーン、UTCからのオフセット'
                          '(時間。JSTなら+9、default: 9.0)')
    ap.add_argument('--site-lat', type=float, default=35.15,
                     help='観測地点の緯度(度、北緯+)。'
                          'default: 35.15(千葉県勝浦市付近)。'
                          '太陽系重心補正への影響は最大でも約21ミリ秒'
                          '(地球半径/光速)なので、今回の観測精度'
                          '(数十秒オーダー)では、多少ずれていても'
                          '実質無視できる。')
    ap.add_argument('--site-lon', type=float, default=140.32,
                     help='観測地点の経度(度、東経+)。default: 140.32')
    ap.add_argument('--site-height', type=float, default=20.0,
                     help='観測地点の標高(m)。default: 20.0')
    args = ap.parse_args()

    jst_dt = datetime.strptime(args.min_time_jst, '%Y-%m-%d %H:%M:%S')
    utc_dt = jst_dt - timedelta(hours=args.tz_offset_hours)

    location = EarthLocation(lat=args.site_lat * u.deg,
                              lon=args.site_lon * u.deg,
                              height=args.site_height * u.m)
    t_utc = Time(utc_dt, scale='utc', location=location)
    coord = SkyCoord(ra=args.ra, dec=args.dec, unit=(u.hourangle, u.deg),
                      frame='icrs')

    # 地心から太陽系重心への光行時間補正(日単位)。
    # 観測地点(地表)の違いによる補正(最大でも地球半径/光速≒21ミリ秒)は、
    # 今回の観測精度(±十数秒程度)に対して無視できるほど小さいため、
    # --site-lat/--site-lon/--site-heightの値が多少ずれていても
    # 結果への影響はほぼ無い。
    ltt_helio = t_utc.light_travel_time(coord, kind='heliocentric')
    hjd_obs = t_utc.jd + ltt_helio.value

    print(f'観測極小時刻(JST): {jst_dt.strftime("%Y-%m-%d %H:%M:%S")}')
    print(f'  → UTC: {utc_dt.strftime("%Y-%m-%d %H:%M:%S")}')
    print(f'  → JD (地心):  {t_utc.jd:.6f}')
    print(f'  → 太陽系重心への光行時間補正: '
          f'{ltt_helio.to(u.minute).value:+.3f} 分')
    print(f'  → HJD (太陽系重心): {hjd_obs:.6f}')
    print()

    cycle = round((hjd_obs - args.epoch) / args.period)
    predicted_hjd = args.epoch + cycle * args.period
    oc_days = hjd_obs - predicted_hjd
    oc_minutes = oc_days * 24 * 60

    print(f'元期(Epoch, HJD): {args.epoch}')
    print(f'周期(Period, days): {args.period}')
    print(f'該当する周回数 E: {cycle}')
    print(f'予報極小時刻(HJD): {predicted_hjd:.6f}')
    print()
    print(f'O-C = {oc_days:+.6f} 日 = {oc_minutes:+.2f} 分')

    if args.min_time_sigma_sec is not None:
        sigma_days = args.min_time_sigma_sec / 86400.0
        sigma_minutes = args.min_time_sigma_sec / 60.0
        print(f'  (極小時刻の1シグマ不確かさ: '
              f'±{sigma_days:.6f}日 = ±{sigma_minutes:.2f}分)')

    return 0


if __name__ == '__main__':
    sys.exit(main())
