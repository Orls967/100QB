# frozen_string_literal: true

module CopulaTailDependence
  module Figures
    module_function

    def write_all(payload, fig_dir, seed:, copula_specs:)
      {
        scatter: write_scatter_fig(fig_dir, seed: seed, copula_specs: copula_specs),
        chi_u: write_chi_fig(payload[:part1_dependence], fig_dir),
        lambda_param: write_lambda_param_fig(fig_dir),
        euler: write_euler_fig(payload[:part2_euler], fig_dir),
        nu_profile: write_nu_profile_fig(payload[:part3_estimation][:misspecification], fig_dir),
        es99_spec: write_es99_spec_fig(payload[:part3_estimation], fig_dir),
        vr_is: write_vr_fig(payload[:part4_importance], fig_dir)
      }
    end

    def write_scatter_fig(fig_dir, seed:, copula_specs:)
      rng = Rng::Sampler.new(seed + 900)
      colors = { gaussian: '#2563eb', t: '#dc2626', clayton: '#059669', gumbel: '#d97706' }
      panels = copula_specs.map do |spec|
        u1 = Array.new(600, 0.0)
        u2 = Array.new(600, 0.0)
        600.times do |i|
          a, b = Copulas.sample_pair(spec[:key], spec[:params], rng)
          u1[i] = a
          u2[i] = b
        end
        { label: spec[:label], u1: u1, u2: u2, color: colors[spec[:key]] }
      end
      out = File.join(fig_dir, 'day005e_copula_scatter.svg')
      Svg.write_scatter_panels(out, title: 'Sampel Kopula Bivariat (N=600 per panel)', panels: panels)
      out
    end

    def write_chi_fig(p1, fig_dir)
      colors = ['#2563eb', '#dc2626', '#059669', '#d97706']
      us = [0.9, 0.99, 0.999]
      series = p1[:copulas].each_with_index.map do |row, idx|
        { label: row[:copula], x: us, y: row[:tail_by_u].map { |r| r[:empirical] }, color: colors[idx], dashed: false }
      end
      out = File.join(fig_dir, 'day005e_tail_chi_u.svg')
      Svg.write_line_chart(out, title: 'Konvergensi Empiris chi(u) vs Ambang u', x_label: 'Ambang u',
                                y_label: 'chi(u)', series: series)
      out
    end

    def write_lambda_param_fig(fig_dir)
      rhos = (1..18).map { |i| i * 0.05 }
      l_g = rhos.map { 0.0 }
      l_t4 = rhos.map { |r| Copulas.exact_upper_tail(:t, rho: r, df: 4.0) }
      l_c = rhos.map { |r| Copulas.exact_lower_tail(:clayton, theta: (2.0 * r) / (1.0 - r + 0.1)) }
      out = File.join(fig_dir, 'day005e_lambda_vs_param.svg')
      Svg.write_line_chart(out, title: 'Koefisien Ketergantungan Ekor vs Parameter', x_label: 'Parameter (rho / skala)',
                                y_label: 'lambda', series: [
                                  { label: 'Gaussian', x: rhos, y: l_g, color: '#2563eb', dashed: true },
                                  { label: 'Student-t (df=4)', x: rhos, y: l_t4, color: '#dc2626', dashed: false },
                                  { label: 'Clayton (lambda_L)', x: rhos, y: l_c, color: '#059669', dashed: false }
                                ])
      out
    end

    def write_euler_fig(p2, fig_dir)
      assets = [1.0, 2.0, 3.0, 4.0, 5.0]
      out = File.join(fig_dir, 'day005e_euler_contribs.svg')
      Svg.write_line_chart(out, title: 'Alokasi Modal Euler ES99 per Aset (d=5)', x_label: 'Indeks Aset',
                                y_label: 'Kontribusi ES99', series: [
                                  { label: 'Gaussian Eksak', x: assets, y: p2[:gaussian_exact][:es_contribs],
                                    color: '#111827', dashed: true },
                                  { label: 'Gaussian MC', x: assets, y: p2[:mc_gaussian][:contribs],
                                    color: '#2563eb', dashed: false },
                                  { label: 'Student-t (df=4)', x: assets, y: p2[:mc_t4][:contribs],
                                    color: '#dc2626', dashed: false },
                                  { label: 'Clayton', x: assets, y: p2[:mc_clayton][:contribs],
                                    color: '#059669', dashed: false }
                                ])
      out
    end

    def write_nu_profile_fig(ms, fig_dir)
      profs = ms[:single_sample][:t_profiles].reject { |r| r[:df].infinite? }
      out = File.join(fig_dir, 'day005e_nu_profile.svg')
      Svg.write_line_chart(out, title: 'Profil Log-Likelihood Kopula Student-t terhadap nu',
                                x_label: 'Derajat Kebebasan nu', y_label: 'Log-Likelihood', series: [
                                  { label: 'Profile LogLik', x: profs.map { |r| r[:df] },
                                    y: profs.map { |r| r[:loglik] }, color: '#2563eb', dashed: false }
                                ])
      out
    end

    def write_es99_spec_fig(p3, fig_dir)
      ms = p3[:misspecification]
      mv = p3[:marginal_vs_copula_es99]
      xs = [1.0, 2.0, 3.0, 4.0]
      ys = [ms[:es99_fitted_gauss_mean], ms[:es99_true_t4], mv[:normal_marginal_t4_copula],
            mv[:t4_marginal_gaussian_copula]]
      out = File.join(fig_dir, 'day005e_es99_misspec.svg')
      Svg.write_line_chart(out, title: 'ES99 Portofolio (1:Gauss, 2:t4, 3:Norm+t4, 4:t4+Gauss)',
                                x_label: 'Spesifikasi Model (1..4)', y_label: 'ES99', series: [
                                  { label: 'ES99 Portofolio', x: xs, y: ys, color: '#dc2626', dashed: false }
                                ])
      out
    end

    def write_vr_fig(p4, fig_dir)
      xs = [1.0, 2.0, 3.0, 4.0]
      ys = p4.map { |r| Math.log10(r[:variance_reduction]) }
      out = File.join(fig_dir, 'day005e_is_variance_reduction.svg')
      Svg.write_line_chart(out, title: 'Log10 Faktor Penghematan Variansi IS (1..4: d2r5, d2r8, d5r5, d5r8)',
                                x_label: 'Skenario', y_label: 'log10(VR Factor)', series: [
                                  { label: 'log10(VR)', x: xs, y: ys, color: '#059669', dashed: false }
                                ])
      out
    end
  end
end
