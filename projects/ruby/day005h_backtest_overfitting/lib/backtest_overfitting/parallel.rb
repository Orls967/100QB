# frozen_string_literal: true

require 'etc'

module BacktestOverfitting
  module Parallel
    CHUNK_SIZE = 500

    module_function

    def default_workers = Etc.nprocessors

    def resolve_workers(workers)
      workers.nil? || workers <= 0 ? default_workers : workers.to_i
    end

    def map_chunks(n_items:, workers: nil, chunk_size: CHUNK_SIZE, &)
      raise ArgumentError, 'n_items harus >= 1' if n_items < 1

      n_chunks = (n_items + chunk_size - 1) / chunk_size
      effective_workers = [resolve_workers(workers), n_chunks].min
      return run_sequential(n_items, n_chunks, chunk_size, &) if effective_workers <= 1 || !Process.respond_to?(:fork)

      run_forked(n_items, n_chunks, chunk_size, effective_workers, &)
    end

    def run_sequential(n_items, n_chunks, chunk_size)
      Array.new(n_chunks) { |chunk_idx| yield(chunk_idx, chunk_length(chunk_idx, n_items, chunk_size)) }
    end
    private_class_method :run_sequential

    def run_forked(n_items, n_chunks, chunk_size, n_workers, &)
      pipes = Array.new(n_workers) { IO.pipe }
      pids = Array.new(n_workers) { |wid| spawn_worker(wid, n_workers, n_items, n_chunks, chunk_size, pipes, &) }
      pipes.each { |(_, writer)| writer.close }
      collect_worker_results(pids, pipes, n_chunks)
    end
    private_class_method :run_forked

    def spawn_worker(wid, n_workers, n_items, n_chunks, chunk_size, pipes, &)
      Process.fork do
        pipes.each_with_index do |(reader, writer), idx|
          reader.close
          writer.close unless idx == wid
        end
        writer = pipes[wid][1]
        writer.binmode
        payload = compute_worker_chunks(wid, n_workers, n_items, n_chunks, chunk_size, &)
        Marshal.dump({ ok: true, results: payload }, writer)
        writer.close
        exit!(0)
      rescue StandardError => e
        Marshal.dump({ ok: false, error: "#{e.class}: #{e.message}" }, writer)
        writer.close
        exit!(1)
      end
    end
    private_class_method :spawn_worker

    def compute_worker_chunks(wid, n_workers, n_items, n_chunks, chunk_size)
      local = []
      c = wid
      while c < n_chunks
        local << [c, yield(c, chunk_length(c, n_items, chunk_size))]
        c += n_workers
      end
      local
    end
    private_class_method :compute_worker_chunks

    def collect_worker_results(pids, pipes, n_chunks)
      results = Array.new(n_chunks)
      errors = []
      pipes.each_with_index do |(reader, _), idx|
        reader.binmode
        raw = reader.read
        reader.close
        _, status = Process.waitpid2(pids[idx])
        record_worker_payload(raw, status, results, errors)
      end
      raise "Worker fork gagal: #{errors.join('; ')}" unless errors.empty?

      results
    end
    private_class_method :collect_worker_results

    def record_worker_payload(raw, status, results, errors)
      return errors << "exit=#{status.exitstatus}" if raw.empty?

      msg = Marshal.load(raw) # rubocop:disable Security/MarshalLoad
      if msg[:ok] && status.success?
        msg[:results].each { |chunk_idx, value| results[chunk_idx] = value }
      else
        errors << (msg[:error] || "exit=#{status.exitstatus}")
      end
    end
    private_class_method :record_worker_payload

    def chunk_length(chunk_idx, n_items, chunk_size)
      [n_items - (chunk_idx * chunk_size), chunk_size].min
    end
    private_class_method :chunk_length
  end
end
