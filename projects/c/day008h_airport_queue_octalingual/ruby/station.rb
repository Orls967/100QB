# frozen_string_literal: true

Server = Struct.new(:status, :enabled, :agent_idx, :extra_sec)

class Station
  attr_accessor :q, :head, :tail, :len, :cap, :has_prio, :c_max, :c_active, :busy_cnt,
                :srv, :arr_cnt, :dep_cnt, :block_cnt, :aban_cnt, :sum_l, :sum_w, :wait_cnt, :busy_slots

  def initialize(c_max, cap, has_prio)
    @q = [0] * MAX_Q
    @head = 0; @tail = 0; @len = 0; @cap = cap; @has_prio = has_prio
    @c_max = c_max; @c_active = c_max; @busy_cnt = 0
    @srv = Array.new(MAX_SRV) { |i| Server.new(0, i < c_max ? 1 : 0, -1, 0) }
    @arr_cnt = 0; @dep_cnt = 0; @block_cnt = 0; @aban_cnt = 0; @sum_l = 0; @sum_w = 0; @wait_cnt = 0; @busy_slots = 0
  end

  def push_tail(idx)
    return if @len >= MAX_Q

    @q[@tail] = idx; @tail = (@tail + 1) % MAX_Q; @len += 1
  end

  def push_head(idx)
    return if @len >= MAX_Q

    @head = (@head - 1 + MAX_Q) % MAX_Q; @q[@head] = idx; @len += 1
  end

  def remove_at(pos)
    idx = @q[(@head + pos) % MAX_Q]
    (pos...(@len - 1)).each { |i| @q[(@head + i) % MAX_Q] = @q[(@head + i + 1) % MAX_Q] }
    @tail = (@tail - 1 + MAX_Q) % MAX_Q; @len -= 1
    idx
  end

  def pop_next(pool, priority_on, prio_br)
    return [-1, prio_br] if @len == 0

    if priority_on != 0 && @has_prio != 0
      @len.times do |i|
        if pool[@q[(@head + i) % MAX_Q]].cls == CLS_BIZ
          return [remove_at(i), prio_br + 1] if i > 0

          break
        end
      end
    end
    res = @q[@head]; @head = (@head + 1) % MAX_Q; @len -= 1
    [res, prio_br]
  end
end
