# frozen_string_literal: true

module CopulaTailDependence
  module Linalg
    module_function

    def cholesky(mat)
      n = mat.length
      l = Array.new(n) { Array.new(n, 0.0) }
      n.times do |i|
        (0..i).each do |j|
          sum = mat[i][j]
          j.times { |k| sum -= l[i][k] * l[j][k] }
          if i == j
            raise ArgumentError, 'Matriks tidak positif tentu' if sum <= 0.0

            l[i][j] = Math.sqrt(sum)
          else
            l[i][j] = sum / l[j][j]
          end
        end
      end
      l
    end

    def dot(a, b)
      sum = 0.0
      i = 0
      n = a.length
      while i < n
        sum += a[i] * b[i]
        i += 1
      end
      sum
    end

    def mat_vec(mat, v)
      mat.map { |row| dot(row, v) }
    end
  end
end
