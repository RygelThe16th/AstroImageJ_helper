#!/usr/bin/env ruby
# frozen_string_literal: true

# tbl_to_csv.rb
#
# AstroImageJ の Measurements Table (.tbl, タブ区切り) を
# プロット用のCSVに変換するスクリプト。
#
# 使い方:
#   ruby tbl_to_csv.rb input.tbl output.csv [options]
#
# 例:
#   # ユリウス日のまま出力
#   ruby tbl_to_csv.rb RZ_CA_Table.tbl rzcas.csv
#
#   # X軸をJST(日本時間)の日時にする
#   ruby tbl_to_csv.rb RZ_CA_Table.tbl rzcas_jst.csv --jst
#
#   # 別の等級列・誤差列を使う場合
#   ruby tbl_to_csv.rb Table2.tbl out.csv --mag-col Source_AMag_T1 --err-col Source_AMag_Err_T1
#
#   # 飽和しているフレームを除外する場合 (Saturated列が0以外の行をスキップ)
#   ruby tbl_to_csv.rb Table.tbl out.csv --exclude-saturated
#
#   # AIJが「No signal for centroid」等でエラーを出した特定フレームを
#   # ファイル名で指定して除外する場合(複数指定可)
#   ruby tbl_to_csv.rb Table.tbl out.csv \
#     --exclude-frame "20241106-214856" --exclude-frame "20241106-215619"
#
#   # 除外したいフレーム名が多いときは、1行1ファイル名のテキストファイルで
#   ruby tbl_to_csv.rb Table.tbl out.csv --exclude-frame-file bad_frames.txt
#
#   # 物理的にあり得ない等級(AAVSOチャートのMax/Min magの範囲外)を除外する場合
#   ruby tbl_to_csv.rb Table.tbl out.csv --min-mag 6.18 --max-mag 7.72
#
# 出力CSVの列:
#   JD (または datetime_jst)  ... X軸
#   mag                        ... 生の等級値
#   plot_y                     ... -mag （等級は小さいほど明るいので、
#                                  そのまま素直にプロットすれば
#                                  明るい方が上に来るようにした値）
#   mag_err                    ... --err-col を指定した場合のみ

require 'optparse'
require 'time'
require 'csv'

# AIJのログの1行から、除外パターン(フレームのファイル名)を取り出す。
# 例: "***ERROR: No signal for centroid in aperture C2 of image
#      Light_rz cas_10.0s_IRCUT_20241106-214856.fit. Multi-Aperture aborted."
# のような行から "Light_rz cas_10.0s_IRCUT_20241106-214856.fit" だけを
# 抜き出す(ファイル名自体に空白やピリオドが含まれるため、単純な空白区切り
# では取れない。".fit"/".fits"の直後にある、文末のピリオドと区別する)。
#
# "of image ... " という形式にマッチしない行(ログの継続行や、無関係な行)
# の場合はnilを返す。ただし、そもそもログ形式ではなく、素のフレーム名や
# 日時の断片が1行1件書かれた従来形式のファイルとの後方互換のため、
# "***"で始まらない行はそのまま(丸ごと)パターンとして扱う。
def extract_exclude_pattern_from_line(line)
  m = line.match(/of image\s+(.+?\.fits?)\.?\s*(?:\s\S.*)?$/i)
  return m[1] if m

  return nil if line.start_with?('***')

  line
end

# 文字列をFloatに変換する。ただしRubyの String#to_f は "NaN" のような
# 数値以外の文字列に対して例外を出さず、黙って 0.0 を返してしまう
# (つまり "NaN".to_f.nan? は false になる)。AIJの.tblには、センタリング
# 失敗などで実際にNaN/Infが書き込まれることがあり、これがそのまま
# 0.0等級として紛れ込む(=測定失敗が実測値であるかのように混入する)
# 実害を確認したため、文字列の段階でNaN/Infを検出してnilを返す。
def safe_parse_float(str)
  return nil if str.nil?

  s = str.strip
  return nil if s.empty?
  return nil if s =~ /\A[+-]?(nan|inf(inity)?)\z/i

  Float(s)
rescue ArgumentError
  nil
end

options = {
  jst: false,
  mag_col: 'Source_AMag_T1',
  jd_col: nil,
  jd_offset: nil,
  jd_col_explicit: false,
  jd_offset_explicit: false,
  err_col: nil,
  exclude_saturated: false,
  sat_col: 'Saturated',
  exclude_frames: [],
  label_col: 'Label',
  min_mag: nil,
  max_mag: nil
}

# --jd-colを指定しなかった場合に自動で探す候補。AIJの'Add astronomical
# data to table'の設定状況次第で、テーブルに実際に出てくるユリウス日系の
# 列名がその都度変わる(JD_UTCが無く、Multi-Aperture素の'J.D.-2400000'
# だけが残っていることがある)ため、よくある候補を順番に試す。
# [列名, その列に対応するオフセット]。
JD_COLUMN_CANDIDATES = [
  ['JD_UTC', 0.0],
  ['HJD_UTC', 0.0],
  ['BJD_TDB', 0.0],
  ['J.D.-2400000', 2400000.0],
  ['JD-2400000', 2400000.0],
  ['J.D.', 0.0]
].freeze

OptionParser.new do |opts|
  opts.banner = 'Usage: ruby tbl_to_csv.rb input.tbl output.csv [options]'

  opts.on('--jst', 'X軸をユリウス日ではなくJST(日本時間)の日時にする') do
    options[:jst] = true
  end

  opts.on('--mag-col NAME', "等級として扱う列名 (default: #{options[:mag_col]})") do |v|
    options[:mag_col] = v
  end

  opts.on('--jd-col NAME',
          'ユリウス日として扱う列名。省略時はJD_UTC/HJD_UTC/BJD_TDB/'\
          '"J.D.-2400000"等、よくある列名を自動的に探す。') do |v|
    options[:jd_col] = v
    options[:jd_col_explicit] = true
  end

  opts.on('--jd-offset N', Float,
          'jd-col オプションで指定した列の値に足すオフセット。' \
          '「J.D.-2400000」のような列を使う場合は 2400000 を指定'\
          '(--jd-colを省略して自動検出させた場合は、検出した列名に'\
          '応じて自動的に設定されるので、通常は指定不要)') do |v|
    options[:jd_offset] = v
    options[:jd_offset_explicit] = true
  end

  opts.on('--err-col NAME', '誤差の列名（指定した場合のみ出力に追加）') do |v|
    options[:err_col] = v
  end

  opts.on('--exclude-saturated',
          "飽和(判定列が0以外)の行を除外する (default判定列: #{options[:sat_col]}、--sat-colで変更可)") do
    options[:exclude_saturated] = true
  end

  opts.on('--sat-col NAME', "飽和判定に使う列名 (default: #{options[:sat_col]})") do |v|
    options[:sat_col] = v
  end

  opts.on('--exclude-frame NAME',
          'AIJのエラーログの行、またはフレーム名(の一部)を指定して、'\
          'そのフレームの行を除外する(部分一致、複数回指定可)。'\
          'ログの行をそのまま渡した場合は自動でファイル名部分だけを'\
          '抜き出す。例: --exclude-frame "20241106-214856"') do |v|
    pattern = extract_exclude_pattern_from_line(v.strip)
    options[:exclude_frames] << pattern if pattern
  end

  opts.on('--exclude-frame-file PATH',
          '除外したいフレーム名を書いたテキストファイルを指定する'\
          '(空行・#始まりの行は無視)。AIJの生のエラーログ(Log.txt等)を'\
          'そのまま渡してもよい('\
          '"of image <ファイル名>." という行から自動的にファイル名だけを'\
          '抜き出す。それ以外の行で***から始まるものは無視し、'\
          'それ以外はファイル名そのものとして扱う、という従来形式との'\
          '後方互換あり)。--exclude-frameと併用可能。') do |v|
    unless File.exist?(v)
      warn "エラー: --exclude-frame-file で指定されたファイルが見つかりません: #{v}"
      exit 1
    end
    File.readlines(v, encoding: 'UTF-8').each do |line|
      line = line.strip
      next if line.empty? || line.start_with?('#')

      pattern = extract_exclude_pattern_from_line(line)
      options[:exclude_frames] << pattern if pattern
    end
  end

  opts.on('--label-col NAME',
          'フレーム名(ファイル名)が入っている列名。--exclude-frame/'\
          '--exclude-frame-fileと組み合わせて使う (default: '\
          "#{options[:label_col]})") do |v|
    options[:label_col] = v
  end

  opts.on('--min-mag MAG', Float,
          '物理的にあり得る最も明るい等級(数値が小さいほど明るい)。'\
          'これより明るい(数値が小さい)行は測定エラーとみなして除外する。'\
          '例: AAVSOチャートのMax magが6.18なら --min-mag 6.18'\
          '（測光誤差の余裕を見て気持ち小さめの値にしても良い）') do |v|
    options[:min_mag] = v
  end

  opts.on('--max-mag MAG', Float,
          '物理的にあり得る最も暗い等級。これより暗い(数値が大きい)行は'\
          '測定エラーとみなして除外する。'\
          '例: AAVSOチャートのMin magが7.72なら --max-mag 7.72') do |v|
    options[:max_mag] = v
  end

  opts.on('-h', '--help', 'このヘルプを表示') do
    puts opts
    exit
  end
end.parse!(ARGV)

input_path, output_path = ARGV

if input_path.nil? || output_path.nil?
  warn 'Usage: ruby tbl_to_csv.rb input.tbl output.csv [options]'
  warn '詳細は --help を参照してください。'
  exit 1
end

unless File.exist?(input_path)
  warn "エラー: 入力ファイルが見つかりません: #{input_path}"
  exit 1
end

# JD (ユリウス日) -> UTC の Time オブジェクト
# JD 2440587.5 = 1970-01-01 00:00:00 UTC (Unixエポック)
def jd_to_utc_time(jd)
  unix_seconds = (jd - 2_440_587.5) * 86_400.0
  Time.at(unix_seconds).utc
end

JST_OFFSET_SEC = 9 * 3600

table = CSV.read(input_path, col_sep: "\t", headers: true)

if options[:jd_col].nil?
  found = JD_COLUMN_CANDIDATES.find { |name, _offset| table.headers.include?(name) }
  if found
    options[:jd_col], candidate_offset = found
    options[:jd_offset] = candidate_offset unless options[:jd_offset_explicit]
  else
    warn 'エラー: ユリウス日の列が見つかりません。次の候補を探しましたが、'\
         'どれもテーブルにありませんでした:'
    warn "  #{JD_COLUMN_CANDIDATES.map(&:first).join(', ')}"
    warn '--jd-col で実際の列名を指定してください。'
    warn "利用可能な列名:\n  #{table.headers.join("\n  ")}"
    exit 1
  end
end
options[:jd_offset] ||= 0.0

[options[:mag_col], options[:jd_col]].each do |col|
  next if table.headers.include?(col)

  warn "エラー: 列 '#{col}' が見つかりません。"
  warn "利用可能な列名:\n  #{table.headers.join("\n  ")}"
  exit 1
end

if options[:err_col] && !table.headers.include?(options[:err_col])
  warn "エラー: 誤差列 '#{options[:err_col]}' が見つかりません。"
  exit 1
end

if options[:exclude_saturated] && !table.headers.include?(options[:sat_col])
  warn "エラー: 飽和判定列 '#{options[:sat_col]}' が見つかりません。"
  warn "利用可能な列名:\n  #{table.headers.join("\n  ")}"
  exit 1
end

if !options[:exclude_frames].empty? && !table.headers.include?(options[:label_col])
  warn "エラー: フレーム名の列 '#{options[:label_col]}' が見つかりません。"
  warn '--label-col で実際の列名を指定してください。'
  warn "利用可能な列名:\n  #{table.headers.join("\n  ")}"
  exit 1
end

written = 0
skipped_saturated = 0
skipped_excluded_frame = 0
skipped_out_of_range = 0
skipped_invalid_value = 0

CSV.open(output_path, 'w') do |out|
  header = [options[:jst] ? 'datetime_jst' : 'JD', 'mag', 'plot_y']
  header << 'mag_err' if options[:err_col]
  out << header

  table.each do |row|
    jd_raw = row[options[:jd_col]]
    mag_raw = row[options[:mag_col]]
    next if jd_raw.nil? || jd_raw.strip.empty?
    next if mag_raw.nil? || mag_raw.strip.empty?

    jd = safe_parse_float(jd_raw)
    mag = safe_parse_float(mag_raw)
    if jd.nil? || mag.nil?
      skipped_invalid_value += 1
      next
    end
    jd += options[:jd_offset]

    if options[:exclude_saturated]
      sat_raw = row[options[:sat_col]]
      sat_value = sat_raw.nil? || sat_raw.strip.empty? ? 0 : sat_raw.to_f
      if sat_value != 0
        skipped_saturated += 1
        next
      end
    end

    unless options[:exclude_frames].empty?
      label_raw = row[options[:label_col]].to_s
      if options[:exclude_frames].any? { |pat| label_raw.include?(pat) }
        skipped_excluded_frame += 1
        next
      end
    end

    if (options[:min_mag] && mag < options[:min_mag]) ||
       (options[:max_mag] && mag > options[:max_mag])
      skipped_out_of_range += 1
      next
    end

    x_value =
      if options[:jst]
        jst_time = jd_to_utc_time(jd) + JST_OFFSET_SEC
        jst_time.strftime('%Y-%m-%d %H:%M:%S')
      else
        format('%.6f', jd)
      end

    line = [x_value, format('%.4f', mag), format('%.4f', -mag)]

    if options[:err_col]
      err_raw = row[options[:err_col]]
      line << (err_raw && !err_raw.strip.empty? ? format('%.4f', err_raw.to_f) : '')
    end

    out << line
    written += 1
  end
end

puts "書き出し完了: #{output_path} (#{written} 行)"
puts "※ NaN/Inf等、数値として不正な行を除外した行: #{skipped_invalid_value} 行" if skipped_invalid_value.positive?
puts "※ 飽和により除外した行: #{skipped_saturated} 行" if options[:exclude_saturated]
puts "※ 指定フレーム名により除外した行: #{skipped_excluded_frame} 行" unless options[:exclude_frames].empty?
puts "※ 物理的にあり得ない等級により除外した行: #{skipped_out_of_range} 行" if options[:min_mag] || options[:max_mag]
puts '※ plot_y 列は -mag です。素直にプロットすれば明るい方が上に来ます。' if written.positive?
