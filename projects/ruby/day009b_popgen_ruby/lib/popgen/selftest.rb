# frozen_string_literal: true

require_relative 'rng'
require_relative 'binomial'
require_relative 'markov'
require_relative 'coalescent'
require_relative 'stats'

module PopGen
  # Single built-in self-test suite (rake selftest).
  module SelfTest
    module_function

    def run
      failures = 0
      puts '[selftest] Starting PopGen Ruby self-test suite...'

      # 1. Binomial sampler on hand-computed small PMF: Bin(n=4, p=0.3)
      # Hand PMF: [0.2401, 0.4116, 0.2646, 0.0756, 0.0081]
      rng = RNG.for_stage(0, 1)
      hand_pmf = [0.2401, 0.4116, 0.2646, 0.0756, 0.0081]
      n_draws = 10_000
      obs = Array.new(5, 0)
      n_draws.times { obs[Binomial.sample(4, 0.3, rng)] += 1 }
      exp = hand_pmf.map { |pr| pr * n_draws }
      gof = Stats.chi2_gof_pooled(obs, exp)
      pmf_err = (0..4).map { |k| (Binomial.pmf(k, 4, 0.3) - hand_pmf[k]).abs }.max
      if pmf_err > 1e-14 || gof[:p_value] <= 0.001
        puts format('  [FAIL] 1. Binomial hand PMF: err=%.3e chi2_p=%.4f', pmf_err, gof[:p_value])
        failures += 1
      else
        puts format('  [OK] 1. Binomial hand PMF Bin(4,0.3): max_err=%.2e, chi2_p=%.4f', pmf_err, gof[:p_value])
      end

      # 2. 2-allele N=2 Wright-Fisher Markov chain computed by hand with Rational:
      # States {0, 1, 2}. From state i=1 (p=1/2): P[1][0]=1/4, P[1][1]=1/2, P[1][2]=1/4.
      # Hand fixation prob u_1 = 1/2; mean absorption time t_1 = 1 / (1 - 1/2) = 2;
      # conditional fixation time w_1 = u_1 / (1 - 1/2) = 1 => t_fix(1) = w_1 / u_1 = 2.
      p_rat = Markov.wf_transition_matrix_rational(2)
      u_rat = Markov.fixation_probs(p_rat, rational: true)
      t_abs = Markov.mean_absorption_times(p_rat, rational: true)
      t_fix = Markov.mean_conditional_fixation_times(p_rat, rational: true)
      ok_wf2 = (p_rat[1] == [Rational(1, 4), Rational(1, 2), Rational(1, 4)]) &&
               (u_rat[1] == Rational(1, 2)) &&
               (t_abs[1] == Rational(2, 1)) &&
               (t_fix[1] == Rational(2, 1))
      if ok_wf2
        puts '  [OK] 2. Hand-computed N=2 WF Rational Markov chain (P, u_1=1/2, t_abs=2, t_fix=2)'
      else
        puts "  [FAIL] 2. N=2 WF Rational chain mismatch: P[1]=#{p_rat[1]} u=#{u_rat} t_abs=#{t_abs}"
        failures += 1
      end

      # 3. Kingman coalescent n=2: E[T_MRCA] = 1 within 3 SE on small replications (R=2000)
      rng3 = RNG.for_stage(0, 3)
      ts = Array.new(2000) { Coalescent.simulate_tree(2, 2.0, rng3)[:t_mrca] }
      m_t = Stats.mean(ts)
      se_t = Stats.se_iid(ts)
      z_t = (m_t - 1.0).abs / se_t
      if z_t > 3.0
        puts format('  [FAIL] 3. Kingman n=2 E[T_MRCA]=%.4f +- %.4f (z=%.2f > 3)', m_t, se_t, z_t)
        failures += 1
      else
        puts format('  [OK] 3. Kingman n=2 E[T_MRCA]=%.4f +- %.4f (z=%.2f <= 3 SE)', m_t, se_t, z_t)
      end

      # 4. Exact Ewens partition sum == 1 for n=6 across theta in {0.5, 1.5, 4.0}
      parts6 = Coalescent.all_partitions(6)
      worst_sum_err = 0.0
      [0.5, 1.5, 4.0].each do |th|
        tot_p = parts6.sum(0.0) { |pt| Coalescent.ewens_partition_prob(pt, 6, th) }
        err = (tot_p - 1.0).abs
        worst_sum_err = err if err > worst_sum_err
      end
      if parts6.length != 11 || worst_sum_err > 1e-14
        puts format('  [FAIL] 4. Ewens n=6 partitions=%d sum_err=%.3e', parts6.length, worst_sum_err)
        failures += 1
      else
        puts format('  [OK] 4. Ewens n=6 partition sum == 1 (11 partitions, max_err=%.2e)', worst_sum_err)
      end

      puts "[selftest] Finished with #{failures} failure(s)."
      failures.zero?
    end
  end
end
