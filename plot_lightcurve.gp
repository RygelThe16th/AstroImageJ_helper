# plot_lightcurve.gp
#
# tbl_to_csv.rb --jst の出力(datetime_jst,mag,plot_y の3列CSV)を
# gnuplotで光度曲線としてプロットするスクリプト。
#
# 使い方:
#   gnuplot -persist -e "csvfile='rzcas_jst.csv'" plot_lightcurve.gp
#
# デフォルトのファイル名を使う場合は csvfile の指定を省略してもよい
# (下の if(!exists...) でデフォルト値を設定している)。
#
# 表示範囲を途中で打ち切りたい場合は、xmin/xmaxを指定する
# (CSVの日付部分も含めたフルの日時文字列で指定すること。
#  -e側でset xrangeを直接指定すると、このファイル内でset xdata time/
#  set timefmtが設定される"前"に評価されてしまいエラーになるため、
#  xmin/xmaxという変数経由でこのファイルの中(timefmt設定より後)で
#  xrangeを適用する仕組みにしている)。
#
#   gnuplot -persist -e "csvfile='rzcas.csv'; xmax='2024-12-06 21:30:00'" \
#       plot_lightcurve.gp

if (!exists("csvfile")) csvfile = 'rzcas_jst.csv'
if (!exists("xmin")) xmin = ''
if (!exists("xmax")) xmax = ''

set datafile separator ','
set key autotitle columnhead

set xdata time
set timefmt '%Y-%m-%d %H:%M:%S'
set format x '%H:%M'
set xlabel 'JST'
set xtics rotate by -45

# xmin/xmaxはここ(timefmt設定より後)で初めて時刻として解釈できる。
if (xmin ne '' && xmax ne '') {
    set xrange [xmin:xmax]
} else { if (xmin ne '') {
    set xrange [xmin:*]
} else { if (xmax ne '') {
    set xrange [*:xmax]
} } }

set ylabel '見かけ等級 (V)'
# 等級は小さいほど明るいので、Y軸そのものを反転させる
set yrange [] reverse

set title 'RZ Cas Light Curve'
set grid

set style data points
set pointsize 0.6

plot csvfile using 1:2 with points pt 7 lc rgb '#1f4e99' title 'V mag'
