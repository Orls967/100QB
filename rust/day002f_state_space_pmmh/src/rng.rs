use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use rand_distr::{Distribution, StandardNormal};

const GOLDEN_GAMMA: u64 = 0x9E37_79B9_7F4A_7C15;
const MIX_MUL1: u64 = 0xBF58_476D_1CE4_E5B9;
const MIX_MUL2: u64 = 0x94D0_49BB_1331_11EB;

pub fn derive_seed(master_seed: u64, stream_id: u64) -> u64 {
    let mut z = master_seed.wrapping_add(stream_id.wrapping_add(1).wrapping_mul(GOLDEN_GAMMA));
    z = (z ^ (z >> 30)).wrapping_mul(MIX_MUL1);
    z = (z ^ (z >> 27)).wrapping_mul(MIX_MUL2);
    z ^ (z >> 31)
}

pub fn make_rng(seed: u64) -> ChaCha8Rng {
    ChaCha8Rng::seed_from_u64(seed)
}

#[inline]
pub fn sample_normal(rng: &mut ChaCha8Rng) -> f64 {
    StandardNormal.sample(rng)
}

#[inline]
pub fn sample_uniform(rng: &mut ChaCha8Rng) -> f64 {
    rng.gen::<f64>()
}

/// Dense lower-triangular Cholesky decomposition `A = L * L^T` for small SPD matrices.
pub fn cholesky_lower(mat: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let n = mat.len();
    if n == 0 {
        return Err("empty matrix in cholesky_lower".to_string());
    }
    let mut lower = vec![vec![0.0; n]; n];
    for i in 0..n {
        if mat[i].len() != n {
            return Err("non-square matrix in cholesky_lower".to_string());
        }
        for j in 0..=i {
            let mut sum = mat[i][j];
            for (&lik, &ljk) in lower[i].iter().zip(&lower[j]).take(j) {
                sum -= lik * ljk;
            }
            if i == j {
                if sum <= 0.0 || !sum.is_finite() {
                    return Err(format!("matrix not positive-definite at diag {i}: {sum}"));
                }
                lower[i][j] = sum.sqrt();
            } else {
                lower[i][j] = sum / lower[j][j];
            }
        }
    }
    Ok(lower)
}

/// Solve `L * x = b` where `L` is lower triangular.
pub fn forward_solve(lower: &[Vec<f64>], rhs: &[f64]) -> Vec<f64> {
    let n = rhs.len();
    let mut x = vec![0.0; n];
    for i in 0..n {
        let mut s = rhs[i];
        for j in 0..i {
            s -= lower[i][j] * x[j];
        }
        x[i] = s / lower[i][i];
    }
    x
}

#[inline]
pub fn normal_log_pdf(x: f64, mean: f64, sd: f64) -> f64 {
    let z = (x - mean) / sd;
    -0.5 * (2.0 * std::f64::consts::PI).ln() - sd.ln() - 0.5 * z * z
}
