# frozen_string_literal: true

Agent = Struct.new(:active, :cls, :dest, :flight_hr, :n_bags, :cur_stn, :extra_sec, :is_bag, :arr_slot, :q_enter_slot, :dep_slot)

def pool_alloc(s)
  return -1 if s.free_stk.empty?

  idx = s.free_stk.pop
  s.pool[idx] = Agent.new(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  idx
end

def pool_free(s, idx)
  return if idx < 0 || idx >= MAX_POOL || s.pool[idx].active == 0

  s.pool[idx].active = 0
  s.free_stk.push(idx)
end
