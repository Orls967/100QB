# frozen_string_literal: true

module GarchFamily
  module Optim
    module_function

    def initial_params(model, returns, dist)
      n = returns.length
      mu = returns.sum(0.0) / n
      var = [returns.sum(0.0) { |r| (r - mu)**2 } / n, 1.0e-8].max
      base = case model
             when Garch11 then { mu: mu, omega: 0.08 * var, alpha: 0.07, beta: 0.85 }
             when GjrGarch11 then { mu: mu, omega: 0.08 * var, alpha: 0.04, gamma: 0.06, beta: 0.85 }
             when Egarch11 then { mu: mu, omega: 0.05 * Math.log(var), alpha: 0.14, gamma: -0.06, beta: 0.95 }
             else raise(ArgumentError, "Model tidak dikenal untuk optimasi: #{model.class}")
             end
      dist.is_a?(Innovations::StudentT) ? base.merge(nu: 8.0) : base
    end

    def fit(model, returns, dist: Innovations::Normal.new, init_params: nil, extra_starts: [],
            max_iter: 120, gtol: 1.0e-5)
      starts = [init_params || initial_params(model, returns, dist)] + extra_starts
      best = nil
      starts.each do |p0|
        res = single_start_fit(model, returns, dist, p0, max_iter, gtol)
        best = res if best.nil? || (res[:loglik].finite? && res[:loglik] > best[:loglik])
        break if res[:converged]
      end
      best
    end

    def single_start_fit(model, returns, dist, p0, max_iter, gtol = 1.0e-5)
      theta0 = model.transform(p0, dist)
      inv_n = 1.0 / returns.length
      sd_r = [sample_std(returns), 1.0e-4].max
      z0 = theta0.dup
      z0[0] = 0.0
      obj = ->(z) { scaled_nll_grad(model, z, theta0[0], sd_r, returns, dist, inv_n) }
      opt = bfgs(z0, obj, max_iter: [max_iter, 160].max, gtol: gtol)
      th_opt = unscale_theta(opt[:x], theta0[0], sd_r)
      ll = -opt[:f] * returns.length
      { params: model.untransform(th_opt, dist), theta: th_opt, loglik: ll, nll: -ll,
        avg_grad_inf: opt[:grad_inf], iterations: opt[:iterations],
        grad_inf: opt[:grad_inf] * returns.length, converged: opt[:converged] }
    end
    private_class_method :single_start_fit

    def sample_std(returns)
      m = returns.sum(0.0) / returns.length
      Math.sqrt(returns.sum(0.0) { |r| (r - m)**2 } / returns.length)
    end
    private_class_method :sample_std

    def unscale_theta(z, mu0, sd_r)
      th = z.dup
      th[0] = mu0 + (sd_r * z[0])
      th
    end
    private_class_method :unscale_theta

    def scaled_nll_grad(model, z, mu0, sd_r, returns, dist, inv_n)
      nll, grad = nll_and_theta_grad(model, unscale_theta(z, mu0, sd_r), returns, dist)
      return [nll, grad] if nll >= 1.0e29

      sc = grad.map { |g| g * inv_n }
      sc[0] *= sd_r
      [nll * inv_n, sc]
    end
    private_class_method :scaled_nll_grad

    def nll_and_theta_grad(model, theta, returns, dist)
      p = model.untransform(theta, dist)
      ll = model.loglik(p, returns, dist)
      return [1.0e30, Array.new(theta.length, 0.0)] unless ll.finite?

      sc = model.score(p, returns, dist)
      grad = chain_neg_grad(model, theta, dist, sc)
      [-ll, grad]
    end

    def chain_neg_grad(model, theta, dist, sc)
      keys = model.param_names(dist)
      h = 1.0e-6
      Array.new(theta.length) do |j|
        tp = theta.dup
        tm = theta.dup
        tp[j] += h
        tm[j] -= h
        pp = model.untransform(tp, dist)
        pm = model.untransform(tm, dist)
        -keys.each_with_index.sum(0.0) { |k, i| sc[i] * ((pp[k] - pm[k]) / (2.0 * h)) }
      end
    end
    private_class_method :chain_neg_grad

    def bfgs(x0, obj, max_iter: 160, gtol: 1.0e-6)
      d = x0.length
      x = x0.dup
      f, g = obj.call(x)
      h_inv = eye(d)
      iters = 0
      step_conv = false
      max_iter.times do |step|
        iters = step + 1
        break if g.map(&:abs).max < gtol

        dir = search_dir(h_inv, g, step)
        alpha, x_new, f_new, g_new = line_search(obj, x, f, g, dir)
        break if alpha.nil?

        s = Array.new(d) { |i| x_new[i] - x[i] }
        y = Array.new(d) { |i| g_new[i] - g[i] }
        df = (f - f_new).abs
        x = x_new
        f = f_new
        g = g_new
        step_conv = step.positive? && s.map(&:abs).max < 1.0e-7 && df < 1.0e-11
        break if step_conv

        h_inv = bfgs_update(h_inv, s, y, d, step)
      end
      ginf = g.map(&:abs).max
      { x: x, f: f, grad: g, grad_inf: ginf, iterations: iters,
        converged: f.finite? && (ginf < [gtol * 10.0, 1.0e-3].max || step_conv) }
    end

    def search_dir(h_inv, g, step)
      dir = mat_vec_neg(h_inv, g)
      dir = g.map(&:-@) if Linalg.dot(dir, g) >= -1.0e-12
      return dir unless step.zero?

      mx = dir.map(&:abs).max
      mx > 1.0 ? dir.map { |v| v / mx } : dir
    end
    private_class_method :search_dir

    def line_search(obj, x, f, g, dir)
      slope = Linalg.dot(g, dir)
      alpha = 1.0
      24.times do
        cand = Array.new(x.length) { |i| x[i] + (alpha * dir[i]) }
        fc, gc = obj.call(cand)
        return [alpha, cand, fc, gc] if fc.finite? && fc <= f + (1.0e-4 * alpha * slope)

        alpha *= 0.5
      end
      [nil, x, f, g]
    end
    private_class_method :line_search

    def bfgs_update(h_inv, s, y, d, step)
      ys = Linalg.dot(y, s)
      sn = Linalg.dot(s, s)
      yn = Linalg.dot(y, y)
      return eye(d) if ys <= 1.0e-12 * Math.sqrt(sn * yn) || !ys.finite?

      h_base = step.zero? ? eye(d, ys / yn) : h_inv
      rho = 1.0 / ys
      hy = Array.new(d) { |i| Linalg.dot(h_base[i], y) }
      yhy = Linalg.dot(y, hy)
      fac = (1.0 + (rho * yhy)) * rho
      Array.new(d) do |i|
        Array.new(d) do |j|
          h_base[i][j] - (rho * ((s[i] * hy[j]) + (hy[i] * s[j]))) + (fac * s[i] * s[j])
        end
      end
    end
    private_class_method :bfgs_update

    def mat_vec_neg(mat, vec)
      Array.new(mat.length) { |i| -Linalg.dot(mat[i], vec) }
    end
    private_class_method :mat_vec_neg

    def eye(d, scale = 1.0)
      Array.new(d) { |i| Array.new(d) { |j| i == j ? scale : 0.0 } }
    end
    private_class_method :eye
  end
end
