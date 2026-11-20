# frozen_string_literal: true

require 'fileutils'
require 'json'

module AmericanPut
  class Runner
    def initialize(params: DEFAULT_PARAMS, seed: DEFAULT_SEED, workers: nil,
                   crr_steps: DEFAULT_CRR_STEPS, lsm_paths: DEFAULT_LSM_PATHS,
                   lsm_steps: DEFAULT_LSM_STEPS, psor_nx: DEFAULT_PSOR_NX, psor_nt: DEFAULT_PSOR_NT)
      @params = params
      @seed = seed
      @workers = workers
      @crr_steps = crr_steps
      @lsm_paths = lsm_paths
      @lsm_steps = lsm_steps
      @psor_nx = psor_nx
      @psor_nt = psor_nt
    end

    def execute
      core = run_core_benchmark
      sweeps = run_all_sweeps
      spot_profile, t_spot = timed { run_spot_profile }
      payload = build_payload(core, sweeps, spot_profile, t_spot)
      {
        payload: payload, crr_sweep: sweeps[:crr], psor_sweep: sweeps[:psor],
        lsm_sweep: sweeps[:lsm], spot_profile: spot_profile
      }
    end

    def write_outputs(result, repo_root:)
      out_dir = File.join(repo_root, 'outputs', 'day002k')
      fig_dir = File.join(repo_root, 'outputs', 'figures')
      FileUtils.mkdir_p([out_dir, fig_dir])
      paths = {
        json: File.join(out_dir, 'results.json'),
        svg_convergence: File.join(fig_dir, 'day002k_american-convergence.svg'),
        svg_spot: File.join(fig_dir, 'day002k_american-spot-profile.svg')
      }
      File.write(paths[:json], "#{JSON.pretty_generate(result[:payload])}\n")
      Svg.write_convergence_svg(paths[:svg_convergence], result[:crr_sweep], result[:psor_sweep],
                                result[:lsm_sweep], result[:payload][:crr][:american_price])
      Svg.write_spot_profile_svg(paths[:svg_spot], result[:spot_profile])
      paths
    end

    def self.find_repo_root(start_dir = Dir.pwd)
      dir = File.expand_path(start_dir)
      until File.directory?(File.join(dir, 'src', 'qb100'))
        raise 'Tidak menemukan root repository 100QB (src/qb100)' if dir == File.dirname(dir)

        dir = File.dirname(dir)
      end
      dir
    end

    private

    def run_core_benchmark
      t0 = monotonic_ms
      bs_eu = Stats.bs_put(@params)
      crr_res, t_crr = timed { Crr.price(@params, n_steps: @crr_steps) }
      lsm_opts = { n_paths: @lsm_paths, n_steps: @lsm_steps, seed: @seed, workers: @workers }
      lsm_res, t_lsm = timed { Lsm.price(@params, **lsm_opts) }
      psor_res, t_psor = timed { Psor.price(@params, n_x: @psor_nx, n_t: @psor_nt) }
      {
        bs_eu: bs_eu, crr: crr_res, lsm: lsm_res, psor: psor_res,
        timings: { crr_ms: t_crr.round(3), lsm_ms: t_lsm.round(3), psor_ms: t_psor.round(3),
                   core_total_ms: (monotonic_ms - t0).round(3) }
      }
    end

    def run_all_sweeps
      crr_sw, t_crr = timed { Crr.convergence_sweep(@params, workers: @workers) }
      psor_sw, t_psor = timed { Psor.convergence_sweep(@params, workers: @workers) }
      lsm_sw, t_lsm = timed { Lsm.convergence_sweep(@params, seed: @seed, workers: @workers) }
      { crr: crr_sw, psor: psor_sw, lsm: lsm_sw, t_crr: t_crr, t_psor: t_psor, t_lsm: t_lsm }
    end

    def run_spot_profile
      SPOT_GRID.map do |spot0|
        prm = @params.with_s0(spot0)
        crr_r = Crr.price(prm, n_steps: @crr_steps)
        lsm_r = Lsm.price(prm, n_paths: @lsm_paths, n_steps: @lsm_steps, seed: @seed, workers: @workers)
        psor_r = Psor.price(prm, n_x: @psor_nx, n_t: @psor_nt)
        SpotComparisonPoint.new(
          s0: spot0, intrinsic: [prm.k - spot0, 0.0].max, bs_european_put: Stats.bs_put(prm),
          crr_american: crr_r.american_price, lsm_american: lsm_r.price, lsm_se: lsm_r.se,
          psor_american: psor_r.american_price, early_exercise_premium: crr_r.early_exercise_premium
        )
      end
    end

    def build_payload(core, sweeps, spot_profile, t_spot)
      {
        params: @params.to_h, seed: @seed, bs_european_put: core[:bs_eu],
        crr: core[:crr].to_h, lsm: core[:lsm].to_h, psor: core[:psor].to_h,
        cross_check: build_cross_check(core[:crr], core[:lsm], core[:psor]),
        convergence: { crr: sweeps[:crr].map(&:to_h), psor: sweeps[:psor].map(&:to_h), lsm: sweeps[:lsm].map(&:to_h) },
        spot_profile: spot_profile.map(&:to_h),
        timings_ms: core[:timings].merge(
          crr_sweep_ms: sweeps[:t_crr].round(3), psor_sweep_ms: sweeps[:t_psor].round(3),
          lsm_sweep_ms: sweeps[:t_lsm].round(3), spot_profile_ms: t_spot.round(3)
        )
      }
    end

    def build_cross_check(crr, lsm, psor)
      {
        crr_vs_psor_diff: crr.american_price - psor.american_price,
        lsm_vs_crr_diff: lsm.price - crr.american_price,
        lsm_vs_crr_z: Stats.z_score(lsm.price - crr.american_price, lsm.se),
        lsm_vs_psor_z: Stats.z_score(lsm.price - psor.american_price, lsm.se)
      }
    end

    def monotonic_ms = Process.clock_gettime(Process::CLOCK_MONOTONIC) * 1000.0

    def timed
      t0 = monotonic_ms
      [yield, monotonic_ms - t0]
    end
  end
end
