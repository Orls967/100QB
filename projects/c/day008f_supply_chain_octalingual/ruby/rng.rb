M64 = 0xFFFFFFFFFFFFFFFF

def splitmix64(x)
  x = (x + 0x9E3779B97F4A7C15) & M64
  z = x
  z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & M64
  z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M64
  [(z ^ (z >> 31)) & M64, x]
end

def fnv1a_bytes(bs)
  h = 0xCBF29CE484222325
  bs.each_byte { |b| h = ((h ^ b) * 0x100000001B3) & M64 }
  h
end

def fnv1a_u32pair(a, b)
  fnv1a_bytes([a, b].pack('V2'))
end

def fnv1a_i64(h, v)
  u = v & M64
  8.times { h = ((h ^ (u & 0xFF)) * 0x100000001B3) & M64; u >>= 8 }
  h
end

def rotl64(x, k)
  ((x << k) | (x >> (64 - k))) & M64
end

class Xoshiro256
  def initialize(seed)
    @s = []
    x = seed & M64
    4.times { z, x = splitmix64(x); @s << z }
  end
  def next_u64
    res = (rotl64((@s[0] + @s[3]) & M64, 23) + @s[0]) & M64
    t = (@s[1] << 17) & M64
    @s[2] ^= @s[0]; @s[3] ^= @s[1]; @s[1] ^= @s[2]; @s[0] ^= @s[3]
    @s[2] ^= t; @s[3] = rotl64(@s[3], 45)
    res
  end
  def rand_int(lo, hi)
    return lo if hi <= lo
    lo + ((next_u64 >> 11) % (hi - lo + 1))
  end
end

class MinHeap
  def initialize; @h = []; end
  def empty?; @h.empty?; end
  def peek_time; @h[0][0]; end
  def ev_lt(a, b); a[0] != b[0] ? a[0] < b[0] : a[1] < b[1]; end
  def push(ev)
    @h << ev
    i = @h.length - 1
    while i > 0
      p = (i - 1) >> 1
      break unless ev_lt(ev, @h[p])
      @h[i] = @h[p]; i = p
    end
    @h[i] = ev
  end
  def pop
    top = @h[0]; last = @h.pop
    unless @h.empty?
      i = 0; n = @h.length
      loop do
        l = 2 * i + 1; r = l + 1; best = i
        best = l if l < n && ev_lt(@h[l], last)
        cmp = best == i ? last : @h[best]
        best = r if r < n && ev_lt(@h[r], cmp)
        break if best == i
        @h[i] = @h[best]; i = best
      end
      @h[i] = last
    end
    top
  end
end
