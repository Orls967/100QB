# frozen_string_literal: true

MASK64 = 0xffffffffffffffff
FNV_OFFSET = 0xcbf29ce484222325
FNV_PRIME = 0x100000001b3

def splitmix64(x)
  x = (x + 0x9e3779b97f4a7c15) & MASK64
  z = x
  z = ((z ^ (z >> 30)) * 0xbf58476d1ce4e5b9) & MASK64
  z = ((z ^ (z >> 27)) * 0x94d049bb133111eb) & MASK64
  [x, (z ^ (z >> 31)) & MASK64]
end

def fnv1a_bytes(bytes_str)
  h = FNV_OFFSET
  bytes_str.each_byte do |b|
    h = ((h ^ b) * FNV_PRIME) & MASK64
  end
  h
end

def fnv1a_u32pair(a, b)
  fnv1a_bytes([a, b].pack("V2"))
end

def fnv1a_i64(h, v)
  [v & MASK64].pack("Q<").each_byte do |b|
    h = ((h ^ b) * FNV_PRIME) & MASK64
  end
  h
end

class Xoshiro256
  def initialize(sd)
    x = sd & MASK64
    @s = Array.new(4) do
      x, val = splitmix64(x)
      val
    end
  end

  def rotl64(v, k)
    ((v << k) | (v >> (64 - k))) & MASK64
  end

  def next_u64
    res = (rotl64((@s[0] + @s[3]) & MASK64, 23) + @s[0]) & MASK64
    t = (@s[1] << 17) & MASK64
    @s[2] ^= @s[0]
    @s[3] ^= @s[1]
    @s[1] ^= @s[2]
    @s[0] ^= @s[3]
    @s[2] ^= t
    @s[3] = rotl64(@s[3], 45)
    res
  end

  def rand_int(lo, hi)
    return lo if hi <= lo

    lo + ((next_u64 >> 11) % (hi - lo + 1))
  end
end
