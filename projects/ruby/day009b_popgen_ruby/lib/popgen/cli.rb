# frozen_string_literal: true

require 'optparse'
require_relative 'io'
require_relative 'selftest'
require_relative 'stages_a'
require_relative 'stages_b'
require_relative 'stages_c'
require_relative 'stages_d'
require_relative 'stages_e'
require_relative 'stages_f'
require_relative 'stages_g'
require_relative 'stages_h'

module PopGen
  # CLI entrypoint supporting selftest, all, group <A..H>, and workers <1|2|4>.
  module CLI
    module_function

    def run(argv)
      cmd = argv.shift
      case cmd
      when 'selftest'
        ok = SelfTest.run
        exit(ok ? 0 : 1)
      when 'workers'
        w = (argv.shift || '4').to_i
        h = StagesH.deterministic_fork_hash(w)
        puts "workers=#{w} fork_fnv1a=#{h}"
      when 'group'
        grp = (argv.shift || 'A').upcase
        ctx = parse_ctx(argv)
        run_group(grp, ctx)
      when 'all'
        ctx = parse_ctx(argv)
        t_start = Time.now
        puts '[popgen] Running selftest first...'
        ctx.selftest_ok = SelfTest.run
        %w[A B C D E F G].each do |g|
          puts "[popgen] Running Group #{g}..."
          run_group(g, ctx)
        end
        puts '[popgen] Running Group H...'
        StagesH.run(ctx, t_start)
        puts format('[popgen] Completed %d stages in %.2fs | Files=%d (+1 JSON) | Det-CSV-FNV1a=%s | Flags=%d',
                    ctx.stages.length, ctx.total_elapsed_sec, ctx.files.length,
                    ctx.combined_csv_hash_hex, ctx.vacuous_flags.length)
        exit(ctx.selftest_ok ? 0 : 1)
      else
        warn 'Usage: popgen <selftest|all|group <A..H>|workers <1|2|4>> [options]'
        exit(1)
      end
    end

    def parse_ctx(argv)
      opts = { scale: 'full', workers: 4, lines: 0, ruby: RUBY_DESCRIPTION, arch: RbConfig::CONFIG['host_cpu'], cores: 4 }
      OptionParser.new do |o|
        o.on('--scale SCALE') { |v| opts[:scale] = v }
        o.on('--workers N', Integer) { |v| opts[:workers] = v }
        o.on('--lines N', Integer) { |v| opts[:lines] = v }
        o.on('--ruby STR') { |v| opts[:ruby] = v }
        o.on('--arch STR') { |v| opts[:arch] = v }
        o.on('--cores N', Integer) { |v| opts[:cores] = v }
      end.parse!(argv)
      IOAudit::RunContext.new(
        scale: opts[:scale], workers: opts[:workers], ruby_lines: opts[:lines],
        ruby_ver: opts[:ruby], arch: opts[:arch], cores: opts[:cores]
      )
    end

    def run_group(grp, ctx)
      case grp
      when 'A' then StagesA.run(ctx)
      when 'B' then StagesB.run(ctx)
      when 'C' then StagesC.run(ctx)
      when 'D' then StagesD.run(ctx)
      when 'E' then StagesE.run(ctx)
      when 'F' then StagesF.run(ctx)
      when 'G' then StagesG.run(ctx)
      when 'H' then StagesH.run(ctx, Time.now)
      else
        raise "Unknown group: #{grp}"
      end
    end
  end
end
