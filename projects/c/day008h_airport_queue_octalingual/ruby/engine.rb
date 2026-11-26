# frozen_string_literal: true

def sim_run(s)
  total_slots = s.rp.days * (s.rp.timescale_div2 != 0 ? SLOTS_PER_DAY * 2 : SLOTS_PER_DAY)
  total_slots.times { |t| sim_step_slot(s, t) }
  rem_pax = 0
  s.pool.each do |ag|
    if ag.active != 0 && ag.is_bag == 0
      rem_pax += 1
      s.sys_sum_sojourn += (total_slots - ag.arr_slot)
    end
  end
  stn_rem = (0...N_STN).sum { |st| st == STN_BELT ? 0 : (s.stn[st].len + s.stn[st].busy_cnt) }
  s.inv_violations += 1 if s.total_arr != s.total_done + s.total_lost + rem_pax || stn_rem != rem_pax
  s.inv_violations += 1 if s.sys_sum_l != s.sys_sum_sojourn
end
