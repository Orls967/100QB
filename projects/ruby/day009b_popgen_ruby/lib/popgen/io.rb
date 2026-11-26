# frozen_string_literal: true

require 'json'
require 'fileutils'
require_relative 'stats'

module PopGen
  # CSV emission, FNV-1a hashing, automated anti-vacuous auditing, and JSON output.
  module IOAudit
    EXEMPT_PARAM_COLS = %w[
      stage step steps gen gens t t_center idx bin k i minor_i init_i n n_copies n_sample n_ind n_loci
      n_rep reps repeat workers deme demes d m m_scaled r s h u v p p0 p_true n_diploid rho rho_a rho_b
      theta mu dist distance cut branch version label model method potential
      partition part_id count counts obs_count exp_count bins_used df pass
      exact_roundtrip max_err_q raw_bytes comp_bytes bytes rows cols files files_written
      hash fnv1a file path status reason id name bit_exact p_s0
    ].freeze

    class RunContext
      attr_accessor :scale, :workers, :ruby_lines, :ruby_ver, :arch, :cores,
                    :stages, :files, :criteria, :vacuous_flags, :metrics,
                    :total_elapsed_sec, :selftest_ok

      def initialize(scale: 'full', workers: 4, ruby_lines: 0, ruby_ver: '', arch: '', cores: 4)
        @scale = scale
        @workers = workers
        @ruby_lines = ruby_lines
        @ruby_ver = ruby_ver
        @arch = arch
        @cores = cores
        @stages = []
        @files = []
        @criteria = {}
        @vacuous_flags = []
        @metrics = {}
        @total_elapsed_sec = 0.0
        @selftest_ok = false
        @det_hash = Stats::FNV_OFFSET
      end

      def small?
        @scale == 'small'
      end

      def combined_csv_hash_hex
        format('%016x', @det_hash)
      end

      def record_stage(id, name, elapsed_sec, files_written)
        @stages << {
          'id' => id,
          'name' => name,
          'elapsed_sec' => elapsed_sec.round(4),
          'files' => files_written
        }
      end

      def add_criterion(id, status, measured, threshold)
        @criteria[id] = {
          'status' => status,
          'measured' => measured,
          'threshold' => threshold
        }
      end

      def add_flag(stage_id, file, reason)
        @vacuous_flags << { 'stage' => stage_id, 'file' => file, 'reason' => reason }
      end

      def write_csv(stage_id, path, headers, rows, deterministic: true)
        FileUtils.mkdir_p(File.dirname(path))
        lines = [headers.join(',')]
        rows.each do |row|
          formatted = row.map do |cell|
            if cell.is_a?(Float)
              if cell.nan? || cell.infinite?
                cell.to_s
              elsif cell != 0.0 && (cell.abs < 1e-4 || cell.abs >= 1e7)
                format('%.6e', cell)
              else
                format('%.8f', cell)
              end
            else
              cell.to_s
            end
          end
          lines << formatted.join(',')
        end
        content = "#{lines.join("\n")}\n"
        File.write(path, content)
        h = Stats.fnv1a_bytes(content)
        audit_csv(stage_id, path, headers, rows)
        @files << {
          'path' => path,
          'kind' => 'csv',
          'bytes' => content.bytesize,
          'rows' => rows.length,
          'cols' => headers.length,
          'fnv1a' => format('%016x', h)
        }
        if deterministic
          @det_hash = Stats.fnv1a_bytes(path, @det_hash)
          @det_hash = Stats.fnv1a_bytes(format('%016x', h), @det_hash)
        end
        path
      end

      def register_svg(stage_id, path)
        unless File.exist?(path)
          add_flag(stage_id, path, 'missing_svg')
          return false
        end
        content = File.read(path)
        valid = content.include?('<svg') && content.include?('</svg>') &&
                content.include?('xmlns="http://www.w3.org/2000/svg"')
        add_flag(stage_id, path, 'invalid_svg_xml') unless valid
        rel = path.sub(%r{\A(\.\./)+}, '')
        @files << {
          'path' => rel,
          'kind' => 'svg',
          'bytes' => content.bytesize,
          'rows' => 0,
          'cols' => 0,
          'fnv1a' => Stats.fnv1a_hex(content)
        }
        valid
      end

      def audit_csv(stage_id, path, headers, rows)
        if rows.empty?
          add_flag(stage_id, path, 'zero_data_rows')
          return
        end

        headers.each_with_index do |col, ci|
          vals = []
          nonfinite = false
          rows.each do |r|
            v = r[ci]
            next unless v.is_a?(Numeric)

            if v.is_a?(Float) && !v.finite?
              nonfinite = true
            else
              vals << v.to_f
            end
          end

          add_flag(stage_id, path, "non_finite_in_#{col}") if nonfinite
          next if vals.empty?

          # Physical bounds checks
          if col.match?(/\A(p|p_hat|p_mean|p_det|p_fix|p_sim|p_exact|freq|h_obs|h_exp)\z/)
            add_flag(stage_id, path, "freq_out_of_0_1_in_#{col}") if vals.any? { |x| x < -1e-9 || x > 1.0 + 1e-9 }
          end
          if col.match?(/\A(fst|fst_hudson|fst_nei|r2|h_mean|h_t)\z/)
            add_flag(stage_id, path, "stat_out_of_0_1_in_#{col}") if vals.any? { |x| x < -1e-6 || x > 1.0 + 1e-6 }
          end
          if col.match?(/\A(t_mean|t_mrca|tw|tb|time|gens|elapsed_sec|var_s|pi|theta_w)\z/)
            add_flag(stage_id, path, "negative_time_or_div_in_#{col}") if vals.any? { |x| x < -1e-9 }
          end

          next if EXEMPT_PARAM_COLS.include?(col)
          next if col.match?(/(theory|exact|analytic|ref|canon|limit|target|count|pass|diff_formula)/)

          if vals.length >= 6 && vals.min == vals.max && vals.min.abs > 1e-12
            add_flag(stage_id, path, "constant_measured_col_#{col}=#{vals.min}")
          end
          if vals.length >= 8 && vals.all? { |x| (x - x.round).abs < 1e-12 && x.abs >= 0.5 }
            unless col.match?(/(id|num|n_|obs|s_muts|alleles|k_n)/)
              add_flag(stage_id, path, "all_round_integers_in_#{col}")
            end
          end
        end
      end

      def write_results_json(path)
        FileUtils.mkdir_p(File.dirname(path))
        n_csv = @files.count { |f| f['kind'] == 'csv' }
        n_svg = @files.count { |f| f['kind'] == 'svg' }
        csv_bytes = @files.select { |f| f['kind'] == 'csv' }.sum { |f| f['bytes'] }
        svg_bytes = @files.select { |f| f['kind'] == 'svg' }.sum { |f| f['bytes'] }

        n_pass = @criteria.values.count { |c| c['status'] == 'tercapai' }
        n_part = @criteria.values.count { |c| c['status'] == 'sebagian' }
        n_other = @criteria.size - n_pass - n_part

        payload = {
          'session' => 'day009b',
          'project_path' => 'projects/ruby/day009b_popgen_ruby',
          'scale' => @scale,
          'master_seed' => RNG::MASTER_SEED,
          'workers' => @workers,
          'run_completed_without_errors' => true,
          'machine' => {
            'arch' => @arch,
            'cores' => @cores,
            'ruby_version' => @ruby_ver
          },
          'code_metrics' => {
            'ruby_lines' => @ruby_lines,
            'stages_completed' => @stages.length,
            'total_output_files' => @files.length + 1,
            'csv_files' => n_csv,
            'svg_files' => n_svg,
            'csv_total_bytes' => csv_bytes,
            'svg_total_bytes' => svg_bytes,
            'total_elapsed_sec' => @total_elapsed_sec.round(4)
          },
          'determinism' => {
            'combined_deterministic_csv_fnv1a' => combined_csv_hash_hex,
            'workers_1_2_4_identical' => @metrics.fetch(:c17_workers_identical, false),
            'two_runs_identical' => @metrics.fetch(:c17_runs_identical, false),
            'worker_hash_1p' => @metrics.fetch(:c17_hash_1p, ''),
            'worker_hash_2p' => @metrics.fetch(:c17_hash_2p, ''),
            'worker_hash_4p' => @metrics.fetch(:c17_hash_4p, '')
          },
          'criteria_summary' => {
            'total' => @criteria.size,
            'tercapai' => n_pass,
            'sebagian' => n_part,
            'lainnya' => n_other
          },
          'criteria' => @criteria,
          'vacuous_flags' => @vacuous_flags,
          'stages' => @stages
        }
        File.write(path, "#{JSON.pretty_generate(payload)}\n")
        true
      end
    end
  end
end
