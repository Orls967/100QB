# frozen_string_literal: true

module CdoTranche
  module ModelRisk
    module_function

    # Membandingkan fair spread (bps) lintas model: :gaussian, :t4, dan :random_recovery.
    def compare_models(n_names: 125, hazard: 0.015, recovery: 0.40, rho: 0.25,
                       n_paths: 20_000, seed: MAIN_SEED, workers: nil)
      %i[gaussian t4 random_recovery].to_h do |model|
        el_curves = simulate_tranche_el_curves(model, n_names, hazard, recovery, rho, n_paths, seed, workers)
        tranches = Models::TRANCHES.each_with_index.map do |tr, idx|
          Pricing.tranche_pv_from_el(el_curves[idx], tr[:attach], tr[:detach]).merge(name: tr[:name])
        end
        [model, tranches]
      end
    end

    def simulate_tranche_el_curves(model, n_names, hazard, recovery, rho, n_paths, seed, workers)
      n_tr = Models::TRANCHES.length
      n_t = Models::PAYMENT_TIMES.length
      c_ps = Models::PAYMENT_TIMES.map { |t| Special.norm_inv(Models.default_prob(hazard, t)) }
      c_ts = Models::PAYMENT_TIMES.map { |t| Special.t_inv(Models.default_prob(hazard, t), 4.0) }
      sr = Math.sqrt(rho)
      s1r = Math.sqrt(1.0 - rho)

      chunks = Parallel.map_chunks(n_items: n_paths, workers: workers, chunk_size: 2_000) do |cidx, len|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, 500 + model_id(model), cidx))
        acc = Array.new(n_tr) { Array.new(n_t, 0.0) }
        len.times { accumulate_path(model, n_names, recovery, sr, s1r, c_ps, c_ts, rng, acc) }
        acc
      end

      merge_curves(chunks, n_tr, n_t, n_paths)
    end

    def model_id(model) = { gaussian: 1, t4: 2, random_recovery: 3 }.fetch(model)
    private_class_method :model_id

    def accumulate_path(model, n_names, recovery, sr, s1r, c_ps, c_ts, rng, acc)
      m = rng.next_normal
      sc = model == :t4 ? Math.sqrt(4.0 / rng.next_chi2(4.0)) : 1.0
      rec = model == :random_recovery ? (recovery + (0.18 * Math.tanh(m))).clamp(0.05, 0.85) : recovery
      lgd_unit = (1.0 - rec) / n_names

      # Simulasi waktu gagal seragam U_i ~ N(0,1) -> default pada t jika Z_i <= thresh(t)
      z_names = Array.new(n_names) { (sr * m) + (s1r * rng.next_normal) }
      Models::PAYMENT_TIMES.each_index do |ti|
        thresh = model == :t4 ? c_ts[ti] / sc : c_ps[ti]
        cnt = z_names.count { |z| z <= thresh }
        loss = cnt * lgd_unit
        Models::TRANCHES.each_with_index do |tr, tri|
          acc[tri][ti] += (loss - tr[:attach]).clamp(0.0, tr[:detach] - tr[:attach])
        end
      end
    end
    private_class_method :accumulate_path

    def merge_curves(chunks, n_tr, n_t, n_paths)
      Array.new(n_tr) do |tri|
        Array.new(n_t) do |ti|
          chunks.sum(0.0) { |c| c[tri][ti] } / n_paths
        end
      end
    end
    private_class_method :merge_curves
  end
end
