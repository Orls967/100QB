# frozen_string_literal: true

module Longevity
  module Parallel
    module_function

    def default_workers
      require 'etc'
      Etc.nprocessors.clamp(1, 8)
    rescue StandardError
      4
    end

    # Maps indices 0...count deterministically across `workers` forked processes.
    # Each item receives `(idx, item_seed)` so results are strictly invariant to `workers`.
    def map_indices(count, workers:, base_seed:, &block)
      return [] if count <= 0

      n_workers = workers.to_i.clamp(1, count)
      return Array.new(count) { |i| block.call(i, Rng.derive_seed(base_seed, i)) } if n_workers == 1

      pipes = []
      pids = []
      n_workers.times do |wid|
        start_idx = (count * wid) / n_workers
        end_idx = (count * (wid + 1)) / n_workers
        rd, wr = IO.pipe
        pid = Process.fork do
          rd.close
          chunk = (start_idx...end_idx).map do |i|
            block.call(i, Rng.derive_seed(base_seed, i))
          end
          Marshal.dump(chunk, wr)
          wr.close
          exit!(0)
        end
        wr.close
        pipes << rd
        pids << pid
      end

      results = []
      pipes.each do |rd|
        results.concat(Marshal.load(rd)) # rubocop:disable Security/MarshalLoad
        rd.close
      end
      pids.each { |pid| Process.wait(pid) }
      results
    end
  end
end
