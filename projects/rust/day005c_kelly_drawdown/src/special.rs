#![allow(clippy::excessive_precision)]
/// Normal PDF, CDF (high-accuracy rational approximation) and inverse CDF.
pub fn norm_pdf(x: f64) -> f64 {
    const INV_SQRT_2PI: f64 = 0.398_942_280_401_432_7;
    INV_SQRT_2PI * (-0.5 * x * x).exp()
}

/// Hart / Cody high-precision erfc(x / sqrt(2)) -> Phi(x) (~1e-15 accuracy).
pub fn norm_cdf(x: f64) -> f64 {
    0.5 * erfc_cody(-x * std::f64::consts::FRAC_1_SQRT_2)
}

fn erfc_cody(x: f64) -> f64 {
    // Cody (1969) rational Chebyshev approximation for erfc(x).
    let ax = x.abs();
    let res = if ax <= 0.46875 {
        const A: [f64; 5] = [
            3.161_123_743_870_565_6,
            1.138_641_541_510_501_6e2,
            3.774_852_376_853_020_2e2,
            3.209_377_589_138_469_5e3,
            1.857_777_061_846_031_5e-1,
        ];
        const B: [f64; 4] = [
            2.360_129_095_234_412_1e1,
            2.440_246_379_344_441_7e2,
            1.282_616_526_077_372_3e3,
            2.844_236_833_439_170_6e3,
        ];
        let ysq = ax * ax;
        let num = (((A[4] * ysq + A[0]) * ysq + A[1]) * ysq + A[2]) * ysq + A[3];
        let den = (((ysq + B[0]) * ysq + B[1]) * ysq + B[2]) * ysq + B[3];
        let erf_val = ax * num / den;
        1.0 - erf_val
    } else if ax <= 4.0 {
        const C: [f64; 9] = [
            5.641_884_969_886_700_9e-1,
            8.883_149_794_388_375_9,
            6.611_919_063_714_162_9e1,
            2.986_351_381_974_001_3e2,
            8.819_522_212_417_690_9e2,
            1.712_047_612_634_070_6e3,
            2.051_078_377_826_071_5e3,
            1.230_339_354_797_997_2e3,
            2.153_115_354_744_038_5e-8,
        ];
        const D: [f64; 8] = [
            1.574_492_611_070_983_5e1,
            1.176_939_508_913_124_9e2,
            5.371_811_018_620_098_6e2,
            1.621_389_574_566_690_2e3,
            3.290_799_235_733_459_6e3,
            4.362_619_090_143_247_2e3,
            3.439_367_674_143_721_6e3,
            1.230_339_354_803_749_4e3,
        ];
        let mut num = C[8];
        for &c in &C[..8] {
            num = num * ax + c;
        }
        let mut den = ax + D[0];
        for &d in &D[1..] {
            den = den * ax + d;
        }
        (-ax * ax).exp() * num / den
    } else if ax < 26.6 {
        const P: [f64; 6] = [
            3.053_266_349_612_323_4e-1,
            3.603_448_999_498_044_4e-1,
            1.257_817_261_112_292_5e-1,
            1.608_378_514_874_227_7e-2,
            6.587_491_615_298_378_0e-4,
            1.631_538_713_730_209_8e-2,
        ];
        const Q: [f64; 5] = [
            2.568_520_192_289_822_2,
            1.872_952_849_923_460_5,
            5.279_051_029_514_284_1e-1,
            6.051_834_131_244_132_0e-2,
            2.335_204_976_268_691_9e-3,
        ];
        let ysq = 1.0 / (ax * ax);
        let mut num = P[5];
        for &p in &P[..5] {
            num = num * ysq + p;
        }
        let mut den = ysq + Q[0];
        for &q in &Q[1..] {
            den = den * ysq + q;
        }
        const FRAC_1_SQRT_PI: f64 = 0.564_189_583_547_756_3;
        (-ax * ax).exp() * (FRAC_1_SQRT_PI - ysq * num / den) / ax
    } else {
        0.0
    };
    if x < 0.0 {
        2.0 - res
    } else {
        res
    }
}

/// Inverse standard normal CDF via bisection + Newton polish.
pub fn inv_norm_cdf(p: f64) -> f64 {
    if p <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }
    let mut lo = -10.0_f64;
    let mut hi = 10.0_f64;
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if norm_cdf(mid) < p {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let mut x = 0.5 * (lo + hi);
    for _ in 0..3 {
        let err = norm_cdf(x) - p;
        let pdf = norm_pdf(x);
        if pdf > 1e-30 {
            x -= err / pdf;
        }
    }
    x
}
