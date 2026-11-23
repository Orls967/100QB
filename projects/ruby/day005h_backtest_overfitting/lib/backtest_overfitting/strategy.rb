# frozen_string_literal: true

module BacktestOverfitting
  module Strategy
    module_function

    def prefix_sums(arr)
      out = Array.new(arr.length + 1, 0.0)
      acc = 0.0
      arr.each_with_index do |v, i|
        acc += v
        out[i + 1] = acc
      end
      out
    end

    def prefix_log_prices(returns)
      px = Array.new(returns.length + 1, 0.0)
      returns.each_with_index do |r, i|
        px[i + 1] = px[i] + Math.log1p([r, -0.999].max)
      end
      prefix_sums(px)
    end

    def build_library
      smas = sma_configs.map { |f, s| SmaCrossover.new(fast: f, slow: s) }
      moms = (1..100).map { |lb| Momentum.new(lookback: lb) }
      mrevs = mr_windows.flat_map { |w| mr_thresholds.map { |z| MeanReversionZScore.new(window: w, threshold: z) } }
      smas + moms + mrevs
    end

    def sma_configs
      fasts = [2, 3, 4, 5, 7, 10, 12, 15, 20, 25]
      slows = [10, 15, 20, 25, 30, 40, 50, 60, 75, 90, 110, 140, 180]
      fasts.flat_map { |f| slows.select { |s| s > f }.map { |s| [f, s] } }.first(100)
    end

    def mr_windows = [5, 8, 10, 12, 15, 18, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 120]
    def mr_thresholds = [0.25, 0.5, 1.0, 1.5, 2.0]

    class Base
      def signal_at(history)
        return 0.0 if history.empty?

        raise NotImplementedError
      end

      def generate_positions(returns)
        Array.new(returns.length) { |t| signal_at(returns[0...t]) }
      end
    end

    class Momentum < Base
      attr_reader :lookback

      def initialize(lookback:)
        super()
        @lookback = lookback
      end

      def name = "MOM(#{@lookback})"

      def signal_at(history)
        return 0.0 if history.length < @lookback

        Strategy.sign_threshold(history.last(@lookback).sum(0.0))
      end

      def generate_positions(returns)
        ps = Strategy.prefix_sums(returns)
        Array.new(returns.length) { |t| t < @lookback ? 0.0 : Strategy.sign_threshold(ps[t] - ps[t - @lookback]) }
      end
    end

    def sign_threshold(val, tol: 1.0e-14)
      return 0.0 if val.abs <= tol

      val.positive? ? 1.0 : -1.0
    end

    class SmaCrossover < Base
      attr_reader :fast, :slow

      def initialize(fast:, slow:)
        super()
        @fast = fast
        @slow = slow
      end

      def name = "SMA(#{@fast},#{@slow})"

      def signal_at(history)
        return 0.0 if history.length < @slow

        px = [0.0]
        history.each { |r| px << (px.last + Math.log1p([r, -0.999].max)) }
        Strategy.sign_threshold((px.last(@fast).sum(0.0) / @fast) - (px.last(@slow).sum(0.0) / @slow))
      end

      def generate_positions(returns)
        pps = Strategy.prefix_log_prices(returns)
        Array.new(returns.length) do |t|
          if t < @slow
            0.0
          else
            ma_f = (pps[t + 1] - pps[t + 1 - @fast]) / @fast
            ma_s = (pps[t + 1] - pps[t + 1 - @slow]) / @slow
            Strategy.sign_threshold(ma_f - ma_s)
          end
        end
      end
    end

    class MeanReversionZScore < Base
      attr_reader :window, :threshold

      def initialize(window:, threshold:)
        super()
        @window = window
        @threshold = threshold
      end

      def name = "MRZ(#{@window},#{@threshold})"

      def signal_at(history)
        return 0.0 if history.length < @window

        win = history.last(@window)
        m = win.sum(0.0) / @window
        sd = Math.sqrt([win.sum(0.0) { |x| (x - m)**2 } / @window, 0.0].max)
        return 0.0 if sd < 1.0e-12

        z_to_pos((history.last - m) / sd)
      end

      def generate_positions(returns)
        ps = Strategy.prefix_sums(returns)
        ps2 = Strategy.prefix_sums(returns.map { |r| r * r })
        Array.new(returns.length) { |t| mr_pos_at(t, returns, ps, ps2) }
      end

      private

      def z_to_pos(z)
        return -1.0 if z > @threshold
        return 1.0 if z < -@threshold

        0.0
      end

      def mr_pos_at(idx, returns, ps, ps2)
        return 0.0 if idx < @window

        m = (ps[idx] - ps[idx - @window]) / @window
        var = ((ps2[idx] - ps2[idx - @window]) / @window) - (m * m)
        return 0.0 if var < 1.0e-24

        z_to_pos((returns[idx - 1] - m) / Math.sqrt(var))
      end
    end
  end
end
