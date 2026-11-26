# frozen_string_literal: true

BP = 10_000
MAX_DAYS = 200
ST_S = 0
ST_E = 1
ST_IA = 2
ST_IS = 3
ST_H = 4
ST_ICU = 5
ST_R = 6
ST_D = 7

Person = Struct.new(
  :age_grp, :hh_id, :sw_id, :status, :variant,
  :state_day, :dur_days, :inf_day, :vax_doses, :vax_day,
  :compliance, :quar_until, :super_flag, :offspring, :ever_inf, :hosp_rej
)

RunParams = Struct.new(
  :n, :days, :beta_bp, :compliance_ovr, :test_ovr, :vax_day_ovr, :vax_cov_ovr,
  :school_ovr, :lock_ovr, :var_day_ovr, :mask_ovr, :contact_scale_bp, :vax_perfect, :interv_on,
  keyword_init: true
)

def default_rp(n, days, beta_bp, interv_on)
  RunParams.new(
    n: n, days: days, beta_bp: beta_bp,
    compliance_ovr: -1, test_ovr: -1, vax_day_ovr: -1, vax_cov_ovr: -1,
    school_ovr: -1, lock_ovr: -1, var_day_ovr: -1, mask_ovr: -1,
    contact_scale_bp: 10_000, vax_perfect: 0, interv_on: interv_on
  )
end

def load_config(path)
  c = {}
  File.readlines(path, chomp: true).each do |line|
    next unless line.include?("=")

    k, v = line.split("=", 2)
    c[k.strip] = v.strip.to_i
  end
  c
end

class Sim
  attr_accessor :cfg, :rp, :n, :hosp_cap, :icu_cap, :hosp_occ, :icu_occ,
                :pop, :n_hh, :n_sw, :hh_off, :hh_mem, :sw_off, :sw_mem, :seen_day,
                :rng, :cand_tgt, :cand_src, :cand_lay, :cand_var,
                :daily_cnt, :daily_inc_age, :daily_deaths_age, :daily_hosp,
                :daily_tests, :daily_quar, :daily_vax, :daily_var,
                :layer_inf, :age_inf, :age_tot, :branches,
                :inv_violations, :total_trans, :total_reinf, :event_hash, :log_enabled,
                :trans, :progs

  def initialize(cfg, rp, log_en)
    @cfg = cfg.dup
    @rp = rp.dup
    @n = rp.n
    @hosp_cap = [(@n * cfg["hosp_beds_per_10k"]) / BP, 1].max
    @icu_cap = [(@n * cfg["icu_beds_per_10k"]) / BP, 1].max
    @hosp_occ = 0
    @icu_occ = 0
    @rng = Array.new(5) do |i|
      x = (cfg["master_seed"] & MASK64) ^ fnv1a_u32pair(100, i)
      _, sd = splitmix64(x)
      Xoshiro256.new(sd)
    end
    n = @n
    @pop = Array.new(n) { Person.new(0, 0, -1, 0, 0, 0, 0, -1, 0, -1, 0, -1, 0, 0, 0, 0) }
    @seen_day = Array.new(n, -1)
    @age_tot = [0, 0, 0]

    hh_cnt = []
    pid = 0
    hh = 0
    while pid < n
      sz = [@rng[0].rand_int(1, 6), n - pid].min
      hh_cnt << sz
      sz.times do
        p = @pop[pid]
        p.hh_id = hh
        u = @rng[0].rand_int(1, 100)
        ag = u <= 22 ? 0 : (u <= 82 ? 1 : 2)
        p.age_grp = ag
        @age_tot[ag] += 1
        p.compliance = @rng[0].rand_int(4000, 9500)
        p.super_flag = @rng[0].rand_int(1, 10) == 1 ? 1 : 0
        pid += 1
      end
      hh += 1
    end
    @n_hh = hh
    @hh_off = Array.new(hh + 1, 0)
    hh.times { |i| @hh_off[i + 1] = @hh_off[i] + hh_cnt[i] }
    @hh_mem = Array.new(n, 0)
    cur = Array.new(hh, 0)
    n.times do |i|
      h = @pop[i].hh_id
      @hh_mem[@hh_off[h] + cur[h]] = i
      cur[h] += 1
    end

    n_schools = (@age_tot[0] / 25) + 1
    n_works = (@age_tot[1] / 15) + 1
    @n_sw = n_schools + n_works
    sw_cnt = Array.new(@n_sw, 0)
    n.times do |i|
      ag = @pop[i].age_grp
      if ag == 0
        sw = @rng[0].rand_int(0, n_schools - 1)
        @pop[i].sw_id = sw
        sw_cnt[sw] += 1
      elsif ag == 1
        sw = n_schools + @rng[0].rand_int(0, n_works - 1)
        @pop[i].sw_id = sw
        sw_cnt[sw] += 1
      end
    end
    @sw_off = Array.new(@n_sw + 1, 0)
    @n_sw.times { |i| @sw_off[i + 1] = @sw_off[i] + sw_cnt[i] }
    @sw_mem = Array.new(n, 0)
    sw_cur = Array.new(@n_sw, 0)
    n.times do |i|
      sw = @pop[i].sw_id
      next if sw < 0

      @sw_mem[@sw_off[sw] + sw_cur[sw]] = i
      sw_cur[sw] += 1
    end

    @cand_tgt = []
    @cand_src = []
    @cand_lay = []
    @cand_var = []
    @daily_cnt = Array.new(MAX_DAYS) { [0] * 8 }
    @daily_inc_age = Array.new(MAX_DAYS) { [0] * 3 }
    @daily_deaths_age = Array.new(MAX_DAYS) { [0] * 3 }
    @daily_hosp = Array.new(MAX_DAYS) { [0] * 4 }
    @daily_tests = Array.new(MAX_DAYS) { [0] * 2 }
    @daily_quar = Array.new(MAX_DAYS, 0)
    @daily_vax = Array.new(MAX_DAYS) { [0] * 2 }
    @daily_var = Array.new(MAX_DAYS) { [0] * 2 }
    @layer_inf = [0, 0, 0]
    @age_inf = [0, 0, 0]
    @branches = [0] * 21
    @inv_violations = 0
    @total_trans = 0
    @total_reinf = 0
    @event_hash = FNV_OFFSET
    @log_enabled = log_en
    @trans = []
    @progs = []
  end
end
