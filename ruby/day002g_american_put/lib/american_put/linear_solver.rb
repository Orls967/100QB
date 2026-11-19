# frozen_string_literal: true

module AmericanPut
  module LinearSolver
    module_function

    # Menyelesaikan sistem linear Ax = b dengan eliminasi Gauss + partial pivoting.
    def solve(matrix_a, rhs_b)
      dim = rhs_b.length
      raise ArgumentError, 'ukuran matriks tidak cocok' unless matrix_a.length == dim

      aug = Array.new(dim) do |i|
        raise ArgumentError, 'matriks tidak persegi' unless matrix_a[i].length == dim

        matrix_a[i].map(&:to_f) + [rhs_b[i].to_f]
      end

      forward_eliminate!(aug, dim)
      back_substitute(aug, dim)
    end

    # Menyelesaikan sistem tridiagonal M x = rhs dengan algoritma Thomas O(n).
    def solve_tridiagonal(lower, diag, upper, rhs)
      size = rhs.length
      c_prime, d_prime = thomas_forward(lower, diag, upper, rhs, size)
      thomas_backward(c_prime, d_prime, size)
    end

    def thomas_forward(lower, diag, upper, rhs, size)
      c_prime = Array.new(size, 0.0)
      d_prime = Array.new(size, 0.0)
      c_prime[0] = upper[0] / diag[0]
      d_prime[0] = rhs[0] / diag[0]

      i = 1
      while i < size
        denom = diag[i] - (lower[i] * c_prime[i - 1])
        raise ZeroDivisionError, 'pivot tridiagonal nol' if denom.abs < 1e-15

        c_prime[i] = upper[i] / denom if i < size - 1
        d_prime[i] = (rhs[i] - (lower[i] * d_prime[i - 1])) / denom
        i += 1
      end
      [c_prime, d_prime]
    end
    private_class_method :thomas_forward

    def thomas_backward(c_prime, d_prime, size)
      sol = Array.new(size, 0.0)
      sol[size - 1] = d_prime[size - 1]
      i = size - 2
      while i >= 0
        sol[i] = d_prime[i] - (c_prime[i] * sol[i + 1])
        i -= 1
      end
      sol
    end
    private_class_method :thomas_backward

    def forward_eliminate!(aug, dim)
      col = 0
      while col < dim
        pivot_row = find_pivot_row(aug, col, dim)
        aug[col], aug[pivot_row] = aug[pivot_row], aug[col] if pivot_row != col
        pivot = aug[col][col]
        raise ZeroDivisionError, 'matriks singular pada eliminasi Gauss' if pivot.abs < 1e-14

        eliminate_below!(aug, col, dim, pivot)
        col += 1
      end
    end
    private_class_method :forward_eliminate!

    def find_pivot_row(aug, col, dim)
      best_row = col
      best_val = aug[col][col].abs
      row = col + 1
      while row < dim
        val = aug[row][col].abs
        if val > best_val
          best_val = val
          best_row = row
        end
        row += 1
      end
      best_row
    end
    private_class_method :find_pivot_row

    def eliminate_below!(aug, col, dim, pivot)
      row = col + 1
      while row < dim
        factor = aug[row][col] / pivot
        if factor != 0.0
          aug[row][col] = 0.0
          k = col + 1
          while k <= dim
            aug[row][k] -= factor * aug[col][k]
            k += 1
          end
        end
        row += 1
      end
    end
    private_class_method :eliminate_below!

    def back_substitute(aug, dim)
      sol = Array.new(dim, 0.0)
      row = dim - 1
      while row >= 0
        sum = aug[row][dim]
        k = row + 1
        while k < dim
          sum -= aug[row][k] * sol[k]
          k += 1
        end
        sol[row] = sum / aug[row][row]
        row -= 1
      end
      sol
    end
    private_class_method :back_substitute
  end
end
