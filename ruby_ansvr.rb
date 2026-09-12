#!/usr/bin/env ruby
# frozen_string_literal: true

# ruby_ansvr.rb
#
# AstroImageJ (および他のnova.astrometry.net互換クライアント) 向けの、
# ローカルsolve-fieldを使ったansvr代替サーバー。
#
# nova.astrometry.netのWeb APIのうち、AIJが使う最小限のサブセット
# (login -> upload -> submissions -> jobs -> wcs_file) だけを実装している。
#
# 前提:
#   - astrometry.net (solve-field) が apt 等でインストール済みで、
#     必要な視野サイズのインデックスファイルも配置済みであること
#   - gem install sinatra puma
#
# 起動:
#   ruby ruby_ansvr.rb
#   (デフォルトで 0.0.0.0:8080 で待ち受ける。ポートは環境変数で変更可)
#
#   ANSVR_PORT=8081 SOLVE_FIELD_BIN=/usr/bin/solve-field ruby ruby_ansvr.rb
#
# デバッグ出力:
#   デフォルトで有効。標準出力(STDOUT)に [ansvr] タグ付きでイベントを
#   逐次出力する(タイムスタンプ、subid、フレームのファイル名、検出星数、
#   solve-fieldの実行コマンドと所要時間など)。
#   ANSVR_DEBUG=0 で抑制できる(エラー等の最低限のログのみ出力)。
#
# AIJ側の設定:
#   WCS > Plate solve with options > Custom Server に
#   http://127.0.0.1:8080 を指定する。
#
# 動作確認のためのcurl例:
#   curl -X POST http://127.0.0.1:8080/api/login
#   curl -X POST http://127.0.0.1:8080/api/upload/ \
#        -F 'request-json={"scale_est":1.19,"scale_err":15,"scale_units":"arcsecperpix"}' \
#        -F 'file=@sample.fits'
#   curl http://127.0.0.1:8080/api/submissions/1
#   curl http://127.0.0.1:8080/api/jobs/1/
#   curl http://127.0.0.1:8080/wcs_file/1 -o result.wcs
#   curl http://127.0.0.1:8080/debug/jobs/1/meta   # フレーム名・検出星数などを確認
#
# 注意:
#   実際にAIJが送ってくるリクエストの細部(フィールド名やポーリングの
#   タイミングなど)と完全に一致しているかは、この場で検証できていない。
#   AIJのWCSログウィンドウ(Watch log的な機能)で実際のやり取りを見ながら
#   微調整が必要になる可能性が高い。

require 'sinatra'
require 'json'
require 'securerandom'
require 'fileutils'
require 'shellwords'
require 'tmpdir'
require 'time'

set :port, ENV.fetch('ANSVR_PORT', 8080).to_i
set :bind, ENV.fetch('ANSVR_BIND', '127.0.0.1')
set :show_exceptions, false

SOLVE_FIELD_BIN = ENV.fetch('SOLVE_FIELD_BIN', 'solve-field')
ASTROMETRY_CFG = ENV.fetch('ASTROMETRY_CFG', '/etc/astrometry.cfg')
PYTHON_BIN = ENV.fetch('PYTHON_BIN', 'python3')
MAKE_XYLS_SCRIPT = ENV.fetch('MAKE_XYLS_SCRIPT', File.join(__dir__, 'make_xyls.py'))
WORK_DIR = ENV.fetch('ANSVR_WORK_DIR', File.join(Dir.tmpdir, 'ruby_ansvr'))
DEBUG = ENV.fetch('ANSVR_DEBUG', '1') != '0'
FileUtils.mkdir_p(WORK_DIR)

# ---------------------------------------------------------------------------
# ロギング
# ---------------------------------------------------------------------------

LOG_MUTEX = Mutex.new

# subid: nil可。ログに [subid=42] のように差し込む(全体の流れを追いやすくするため)。
def log(msg, subid: nil, level: :info)
  return if level == :debug && !DEBUG

  tag = subid ? "[ansvr][subid=#{subid}]" : '[ansvr]'
  line = "#{Time.now.strftime('%Y-%m-%d %H:%M:%S.%3N')} #{tag} #{msg}"
  LOG_MUTEX.synchronize do
    stream = level == :error ? $stderr : $stdout
    stream.puts(line)
    stream.flush
  end
end

# 起動時にパス解決の状況を確認しておく(__dir__はsymlinkを解決するので、
# ruby_ansvr.rb自体をsymlinkで運用していても、実体があるディレクトリを指す)。
log("起動: port=#{ENV.fetch('ANSVR_PORT', 8080)} bind=#{ENV.fetch('ANSVR_BIND', '127.0.0.1')}")
log("SOLVE_FIELD_BIN=#{SOLVE_FIELD_BIN}")
log("ASTROMETRY_CFG=#{ASTROMETRY_CFG} (exists: #{File.exist?(ASTROMETRY_CFG)})")
log("MAKE_XYLS_SCRIPT=#{MAKE_XYLS_SCRIPT} (exists: #{File.exist?(MAKE_XYLS_SCRIPT)})")
log("WORK_DIR=#{WORK_DIR}")
log("DEBUG=#{DEBUG} (ANSVR_DEBUG=0 で抑制可能)")

unless File.exist?(MAKE_XYLS_SCRIPT)
  log("警告: MAKE_XYLS_SCRIPTが見つかりません。テキストxylist(AIJの検出結果)を" \
      '受け取った場合の変換に失敗します。', level: :error)
end

# job_id => { status:, wcs_path:, log: }
JOBS = {}
# subid => { jobs: [job_id], processing_finished: Float|nil }
SUBMISSIONS = {}
NEXT_ID_MUTEX = Mutex.new
NEXT_ID = { value: 0 }

def next_id
  NEXT_ID_MUTEX.synchronize do
    NEXT_ID[:value] += 1
    NEXT_ID[:value]
  end
end

def looks_like_fits?(path)
  head = File.binread(path, 6)
  head == 'SIMPLE'
rescue StandardError
  false
end

# プレーンテキストのxylist(AIJが送ってくる検出結果)の行数(=星数らしきもの)を数える。
# 空行や、数値2つ以上を持たない行は除外する。
def count_text_stars(path)
  File.readlines(path).count do |line|
    parts = line.strip.split
    parts.length >= 2 && parts[0..1].all? { |p| p =~ /\A-?\d+(\.\d+)?\z/ }
  end
rescue StandardError => e
  log("count_text_stars失敗: #{e.class}: #{e.message}", level: :error)
  -1
end

def build_solve_field_args_for_image(input_path, job_dir, opts)
  args = [
    SOLVE_FIELD_BIN,
    '--overwrite',
    '--no-plots',
    '-D', job_dir
  ]

  if ASTROMETRY_CFG
    args += ['--backend-config', ASTROMETRY_CFG]
  end

  if opts['scale_est'] && opts['scale_err']
    est = opts['scale_est'].to_f
    err_pct = opts['scale_err'].to_f
    low = est * (1 - err_pct / 100.0)
    high = est * (1 + err_pct / 100.0)
    units = opts['scale_units'] || 'arcsecperpix'
    args += ['-u', units, '-L', low.to_s, '-H', high.to_s]
  elsif opts['scale_lower'] && opts['scale_upper']
    units = opts['scale_units'] || 'arcsecperpix'
    args += ['-u', units, '-L', opts['scale_lower'].to_s, '-H', opts['scale_upper'].to_s]
  end

  if opts['center_ra'] && opts['center_dec']
    args += ['--ra', opts['center_ra'].to_s, '--dec', opts['center_dec'].to_s]
    args += ['--radius', (opts['radius'] || 1).to_s]
  end

  args << input_path
  args
end

def build_solve_field_args_for_xylist(xyls_path, job_dir, opts)
  args = [
    SOLVE_FIELD_BIN,
    '--overwrite',
    '--no-plots',
    '-D', job_dir
  ]

  if ASTROMETRY_CFG
    args += ['--backend-config', ASTROMETRY_CFG]
  end

  width = opts['image_width']
  height = opts['image_height']
  raise 'image_width/image_height が request-json にありません' unless width && height

  args += ['-w', width.to_s, '-e', height.to_s, '-X', 'X', '-Y', 'Y']

  if opts['scale_est'] && opts['scale_err']
    est = opts['scale_est'].to_f
    err_pct = opts['scale_err'].to_f
    low = est * (1 - err_pct / 100.0)
    high = est * (1 + err_pct / 100.0)
    units = opts['scale_units'] || 'arcsecperpix'
    args += ['-u', units, '-L', low.to_s, '-H', high.to_s]
  elsif opts['scale_lower'] && opts['scale_upper']
    units = opts['scale_units'] || 'arcsecperpix'
    args += ['-u', units, '-L', opts['scale_lower'].to_s, '-H', opts['scale_upper'].to_s]
  end

  if opts['center_ra'] && opts['center_dec']
    args += ['--ra', opts['center_ra'].to_s, '--dec', opts['center_dec'].to_s]
    args += ['--radius', (opts['radius'] || 1).to_s]
  end

  args << xyls_path
  args
end

def run_solve_field(args, subid: nil)
  cmd = args.map { |a| Shellwords.escape(a.to_s) }.join(' ')
  log("solve-field実行: #{cmd}", subid: subid, level: :debug)
  started = Time.now
  log_output = `#{cmd} 2>&1`
  elapsed = Time.now - started
  log("solve-field終了: #{format('%.2f', elapsed)}秒", subid: subid, level: :debug)
  [cmd, log_output, elapsed]
end

post '/api/login' do
  content_type :json
  log('login')
  { status: 'success', session: SecureRandom.hex(8) }.to_json
end

post %r{/api/upload/?} do
  content_type :json

  raw_json = params['request-json'] || params['request_json']
  opts =
    begin
      raw_json ? JSON.parse(raw_json) : {}
    rescue JSON::ParserError => e
      log("request-jsonのパースに失敗: #{e.message}", level: :error)
      {}
    end

  file_param = params[:file]
  if file_param.nil? || file_param[:tempfile].nil?
    log('アップロードエラー: ファイルが添付されていません', level: :error)
    status 400
    next { status: 'error', errormessage: 'no file uploaded' }.to_json
  end

  subid = next_id
  job_id = subid # このサーバーでは1submission = 1jobとして扱う

  original_filename = file_param[:filename].to_s
  log("アップロード受信: filename=#{original_filename.inspect} opts=#{opts.inspect}", subid: subid)

  job_dir = File.join(WORK_DIR, subid.to_s)
  FileUtils.mkdir_p(job_dir)
  input_path = File.join(job_dir, 'input.fits')
  FileUtils.cp(file_param[:tempfile].path, input_path)

  JOBS[job_id] = { status: 'solving', wcs_path: nil, log: '' }
  SUBMISSIONS[subid] = { jobs: [job_id], processing_finished: nil }

  # デバッグ用: 受け取ったrequest-json・元ファイル名・ファイルの先頭バイトを保存しておく。
  # original_filename.txt があれば、後から「subid Nはどのフレームだったか」を
  # 逆算せずに直接特定できる。
  File.write(File.join(job_dir, 'request.json'), raw_json.to_s)
  File.write(File.join(job_dir, 'original_filename.txt'), original_filename)
  File.write(File.join(job_dir, 'file_head.txt'), File.binread(input_path, 200).inspect)

  is_fits_image = looks_like_fits?(input_path)
  # プレーンテキスト(AIJの検出結果)の場合は、その場で星数を数えてファイルに残す。
  # 「AIJが検出した数」と「サーバーが受け取った数」を突き合わせるための材料。
  unless is_fits_image
    star_count = count_text_stars(input_path)
    File.write(File.join(job_dir, 'star_count.txt'), star_count.to_s)
    log("受信した星数: #{star_count}", subid: subid)
  end

  # ここでは即座にレスポンスを返し、実際の解決処理はバックグラウンドの
  # スレッドで行う(nova.astrometry.net本来の非同期プロトコルに合わせる)。
  # AIJはこの後 /api/submissions/<subid> をポーリングして完了を待つ。
  Thread.new do
    begin
      if is_fits_image
        log('入力はFITS画像と判定', subid: subid, level: :debug)
        args = build_solve_field_args_for_image(input_path, job_dir, opts)
      else
        log('入力はプレーンテキスト(xylist候補)と判定', subid: subid, level: :debug)
        xyls_path = File.join(job_dir, 'input.xyls')
        py_cmd = [PYTHON_BIN, MAKE_XYLS_SCRIPT, input_path, xyls_path]
                  .map { |a| Shellwords.escape(a.to_s) }.join(' ')
        log("make_xyls.py実行: #{py_cmd}", subid: subid, level: :debug)
        py_log = `#{py_cmd} 2>&1`
        File.write(File.join(job_dir, 'make_xyls.log'), py_log)

        if File.exist?(xyls_path)
          args = build_solve_field_args_for_xylist(xyls_path, job_dir, opts)
        else
          log("make_xyls.py失敗:\n#{py_log}", subid: subid, level: :error)
          JOBS[job_id][:status] = 'failure'
          JOBS[job_id][:log] = "make_xyls.py failed:\n#{py_log}"
          SUBMISSIONS[subid][:processing_finished] = Time.now.to_f
          next
        end
      end

      _cmd, solve_log, elapsed = run_solve_field(args, subid: subid)
      JOBS[job_id][:log] = solve_log

      wcs_path = File.join(job_dir, 'input.wcs')
      if File.exist?(wcs_path)
        JOBS[job_id][:status] = 'success'
        JOBS[job_id][:wcs_path] = wcs_path
        log("解決成功 (#{format('%.2f', elapsed)}秒)", subid: subid)
      else
        JOBS[job_id][:status] = 'failure'
        log("解決失敗 (#{format('%.2f', elapsed)}秒)", subid: subid, level: :error)
      end
    rescue StandardError => e
      JOBS[job_id][:status] = 'failure'
      JOBS[job_id][:log] = "internal error: #{e.class}: #{e.message}"
      log("内部エラー: #{e.class}: #{e.message}\n#{e.backtrace&.first(5)&.join("\n")}",
          subid: subid, level: :error)
    ensure
      SUBMISSIONS[subid][:processing_finished] = Time.now.to_f
    end
  end

  { status: 'success', subid: subid }.to_json
end

get %r{/api/submissions/(\d+)} do |subid|
  content_type :json
  sub = SUBMISSIONS[subid.to_i]
  halt 404, { status: 'error' }.to_json unless sub

  {
    jobs: sub[:jobs],
    job_calibrations: sub[:jobs].map { |j| [j, j] },
    processing_finished: sub[:processing_finished] || false
  }.to_json
end

get %r{/api/jobs/(\d+)/?} do |job_id|
  content_type :json
  job = JOBS[job_id.to_i]
  halt 404, { status: 'error' }.to_json unless job

  { status: job[:status] }.to_json
end

get %r{/wcs_file/(\d+)} do |job_id|
  job = JOBS[job_id.to_i]
  halt 404 unless job && job[:wcs_path] && File.exist?(job[:wcs_path])
  send_file job[:wcs_path], type: 'application/octet-stream'
end

# デバッグ用: solve-fieldの生ログを見たいとき
get %r{/debug/jobs/(\d+)/log} do |job_id|
  content_type 'text/plain'
  job = JOBS[job_id.to_i]
  halt 404, 'no such job' unless job
  job[:log]
end

# デバッグ用: そのジョブでAIJから受け取ったrequest-jsonの中身を見たいとき
get %r{/debug/jobs/(\d+)/request} do |job_id|
  content_type 'text/plain'
  path = File.join(WORK_DIR, job_id, 'request.json')
  halt 404, 'no such job' unless File.exist?(path)
  File.read(path)
end

# デバッグ用: そのジョブの元フレーム名・検出星数・ジョブの状態をまとめて見たいとき。
# subidからフレームを逆算する必要がないように、これを最初に見るのがおすすめ。
get %r{/debug/jobs/(\d+)/meta} do |job_id|
  content_type :json
  job_dir = File.join(WORK_DIR, job_id)
  halt 404, { status: 'error', errormessage: 'no such job dir' }.to_json unless Dir.exist?(job_dir)

  original_filename_path = File.join(job_dir, 'original_filename.txt')
  star_count_path = File.join(job_dir, 'star_count.txt')

  {
    job_id: job_id.to_i,
    original_filename: (File.read(original_filename_path) if File.exist?(original_filename_path)),
    star_count: (File.read(star_count_path).to_i if File.exist?(star_count_path)),
    status: JOBS[job_id.to_i] && JOBS[job_id.to_i][:status],
    job_dir: job_dir
  }.to_json
end
