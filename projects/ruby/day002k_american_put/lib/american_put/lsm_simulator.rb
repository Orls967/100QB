# frozen_string_literal: true

module AmericanPut
  module LsmSimulator
    module_function

    def simulate_all_paths(params, n_paths, n_steps, step_dt, seed, antithetic, workers)
      chunks = Parallel.map_chunks(n_items: n_paths, workers: workers, chunk_size: Parallel::CHUNK_SIZE) do |c_idx, len|
        c_seed = Rng.derive_chunk_seed(seed, c_idx)
        simulate_chunk_packed(params, len, n_steps, step_dt, c_seed, antithetic)
      end
      unpack_columns(chunks, n_paths, n_steps)
    end

    def simulate_chunk_packed(params, len, n_steps, step_dt, seed, antithetic)
      rng = Rng::NormalGenerator.new(seed)
      drift = (params.r - (0.5 * params.sigma * params.sigma)) * step_dt
      vol = params.sigma * Math.sqrt(step_dt)
      normals = Array.new(n_steps, 0.0)
      flat = Array.new(len * n_steps, 0.0)

      if antithetic
        fill_antithetic_chunk!(flat, len, n_steps, params.s0, drift, vol, normals, rng)
      else
        fill_standard_chunk!(flat, len, n_steps, params.s0, drift, vol, normals, rng)
      end

      flat.pack('d*')
    end

    def fill_standard_chunk!(flat, len, n_steps, spot0, drift, vol, normals, rng)
      p_idx = 0
      while p_idx < len
        rng.fill_normals!(normals)
        spot = spot0
        base = p_idx * n_steps
        step = 0
        while step < n_steps
          spot *= Math.exp(drift + (vol * normals[step]))
          flat[base + step] = spot
          step += 1
        end
        p_idx += 1
      end
    end

    def fill_antithetic_chunk!(flat, len, n_steps, spot0, drift, vol, normals, rng)
      half = len / 2
      p_idx = 0
      while p_idx < half
        rng.fill_normals!(normals)
        fill_pair_steps!(flat, p_idx * 2 * n_steps, n_steps, spot0, drift, vol, normals)
        p_idx += 1
      end
    end

    def fill_pair_steps!(flat, base, n_steps, spot0, drift, vol, normals)
      s_pos = spot0
      s_neg = spot0
      step = 0
      while step < n_steps
        shock = vol * normals[step]
        s_pos *= Math.exp(drift + shock)
        s_neg *= Math.exp(drift - shock)
        flat[base + step] = s_pos
        flat[base + n_steps + step] = s_neg
        step += 1
      end
    end

    def unpack_columns(packed_chunks, n_paths, n_steps)
      cols = Array.new(n_steps) { Array.new(n_paths, 0.0) }
      offset = 0
      packed_chunks.each do |raw|
        vals = raw.unpack('d*')
        chunk_paths = vals.length / n_steps
        copy_chunk_to_columns!(cols, vals, offset, chunk_paths, n_steps)
        offset += chunk_paths
      end
      cols
    end

    def copy_chunk_to_columns!(cols, vals, offset, chunk_paths, n_steps)
      p_idx = 0
      while p_idx < chunk_paths
        base = p_idx * n_steps
        dst = offset + p_idx
        step = 0
        while step < n_steps
          cols[step][dst] = vals[base + step]
          step += 1
        end
        p_idx += 1
      end
    end
  end
end
