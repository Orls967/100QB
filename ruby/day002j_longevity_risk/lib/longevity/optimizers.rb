# frozen_string_literal: true

module Longevity
  module Optimizers
    module_function

    def nelder_mead(x_init, step:, tol: 1e-12, max_iter: 2000, &objective)
      n = x_init.size
      simplex = init_simplex(x_init, step)
      vals = simplex.map { |pt| objective.call(pt) }

      max_iter.times do
        order = (0..n).sort_by { |i| vals[i] }
        simplex = order.map { |i| simplex[i] }
        vals = order.map { |i| vals[i] }
        break if converged?(simplex, vals, n, tol)

        step_simplex!(simplex, vals, compute_centroid(simplex, n), n, &objective)
      end

      { x: simplex[0], f: vals[0] }
    end

    def init_simplex(x_init, step)
      pts = [x_init.dup]
      x_init.size.times do |i|
        pt = x_init.dup
        pt[i] += step[i]
        pts << pt
      end
      pts
    end

    def converged?(simplex, vals, dim, tol)
      spread = (0..dim).map { |i| (vals[i] - vals[0]).abs }.max
      pspread = (1..dim).map { |i| (0...dim).map { |d| (simplex[i][d] - simplex[0][d]).abs }.max }.max
      spread < tol && pspread < tol
    end

    def compute_centroid(simplex, dim)
      Array.new(dim) do |d|
        s = 0.0
        dim.times { |i| s += simplex[i][d] }
        s / dim.to_f
      end
    end

    def combine(vec_a, vec_b, alpha, beta)
      Array.new(vec_a.size) { |d| (alpha * vec_a[d]) + (beta * vec_b[d]) }
    end

    def step_simplex!(simplex, vals, centroid, dim, &objective)
      xr = combine(centroid, simplex[dim], 2.0, -1.0)
      fr = objective.call(xr)
      return expand_or_reflect!(simplex, vals, centroid, xr, fr, dim, &objective) if fr < vals[0]

      if fr < vals[dim - 1]
        simplex[dim] = xr
        vals[dim] = fr
        return
      end

      contract_or_shrink!(simplex, vals, centroid, xr, fr, dim, &objective)
    end

    def expand_or_reflect!(simplex, vals, centroid, x_ref, f_ref, dim, &objective)
      xe = combine(centroid, x_ref, -1.0, 2.0)
      fe = objective.call(xe)
      simplex[dim] = fe < f_ref ? xe : x_ref
      vals[dim] = [fe, f_ref].min
    end

    def contract_or_shrink!(simplex, vals, centroid, x_ref, f_ref, dim, &objective)
      target = f_ref < vals[dim] ? x_ref : simplex[dim]
      xc = combine(centroid, target, 0.5, 0.5)
      fc = objective.call(xc)
      if fc < [f_ref, vals[dim]].min
        simplex[dim] = xc
        vals[dim] = fc
        return
      end

      (1..dim).each do |i|
        simplex[i] = combine(simplex[0], simplex[i], 0.5, 0.5)
        vals[i] = objective.call(simplex[i])
      end
    end
  end
end
