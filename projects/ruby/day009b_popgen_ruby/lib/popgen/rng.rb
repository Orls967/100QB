# frozen_string_literal: true

module PopGen
  # Deterministic RNG and process-forking helpers derived from master seed 20261127.
  module RNG
    MASTER_SEED = 20_261_127
    MASK64 = 0xFFFFFFFFFFFFFFFF

    module_function

    def splitmix64(state)
      z = (state + 0x9E3779B97F4A7C15) & MASK64
      z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
      z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK64
      (z ^ (z >> 31)) & MASK64
    end

    def derive_seed(stage_id, task_id = 0, master = MASTER_SEED)
      s = (master ^ ((stage_id.to_i + 1) * 0x9E3779B97F4A7C15)) & MASK64
      s = splitmix64(s)
      s = (s ^ ((task_id.to_i + 1) * 0xBF58476D1CE4E5B9)) & MASK64
      splitmix64(s)
    end

    def for_stage(stage_id, task_id = 0, master = MASTER_SEED)
      Random.new(derive_seed(stage_id, task_id, master))
    end

    # Standard normal N(0, 1) via polar Box-Muller on Random#rand
    def rand_normal(rng)
      loop do
        u1 = 2.0 * rng.rand - 1.0
        u2 = 2.0 * rng.rand - 1.0
        s = u1 * u1 + u2 * u2
        return u1 * Math.sqrt(-2.0 * Math.log(s) / s) if s > 1e-300 && s < 1.0
      end
    end

    # Deterministic parallel map using Process.fork + IO.pipe (never Ractor).
    # Merges task outputs in strict task-index order so 1/2/4 processes are bitwise identical.
    def parallel_map(tasks, workers: 1)
      n = tasks.length
      return [] if n.zero?

      w = workers.to_i
      w = 1 if w < 1
      w = n if w > n

      if w == 1
        return tasks.each_with_index.map { |task, idx| yield(task, idx) }
      end

      results = Array.new(n)
      slices = (0...n).each_slice((n.to_f / w).ceil).to_a
      children = []

      slices.each do |idx_slice|
        rd, wr = IO.pipe
        pid = Process.fork do
          rd.close
          local = idx_slice.map { |idx| [idx, yield(tasks[idx], idx)] }
          wr.binmode
          wr.write(Marshal.dump(local))
          wr.close
          exit!(0)
        end
        wr.close
        children << [pid, rd]
      end

      children.each do |pid, rd|
        rd.binmode
        payload = rd.read
        rd.close
        Process.waitpid(pid)
        raise "Forked worker #{pid} failed" unless $?.success?

        Marshal.load(payload).each do |idx, val| # rubocop:disable Security/MarshalLoad
          results[idx] = val
        end
      end

      results
    end
  end
end
