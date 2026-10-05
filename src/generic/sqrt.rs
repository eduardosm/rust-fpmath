use crate::traits::{Float, Int as _};

pub(crate) trait Sqrt: Float {
    fn sqrt_finite(x: Self) -> Self;
}

pub(crate) fn sqrt<F: Sqrt>(x: F) -> F {
    let raw = x.to_raw();
    if raw != F::Raw::ZERO && raw < F::EXP_MASK {
        // x is positive, finite and non-zero
        F::sqrt_finite(x)
    } else if (raw & !F::SIGN_MASK) == F::Raw::ZERO {
        // sqrt(±0) = ±0
        x
    } else if x.is_sign_negative() {
        // x < 0, sqrt(x) = NaN
        F::NAN
    } else {
        // propagate infinity or NaN
        x
    }
}

// GENERATE: rsqrt_table
// RSQRT_TBL[i] ~= 1 / sqrt(m) (0.16 fixed point), minimizing the maximum
// relative error for m in [a, b], with
// * a = 2^(i >> 6) * (1 + (i & 63) / 64)
// * b = 2^(i >> 6) * (1 + ((i & 63) + 1) / 64)
static RSQRT_TBL: [u16; 128] = [
    0xFF02, 0xFD0E, 0xFB25, 0xF947, 0xF773, 0xF5AA, 0xF3EA, 0xF234, 0xF087, 0xEEE3, 0xED47, 0xEBB3,
    0xEA27, 0xE8A3, 0xE727, 0xE5B2, 0xE443, 0xE2DC, 0xE17A, 0xE020, 0xDECB, 0xDD7D, 0xDC34, 0xDAF1,
    0xD9B3, 0xD87B, 0xD748, 0xD61A, 0xD4F1, 0xD3CD, 0xD2AD, 0xD192, 0xD07B, 0xCF69, 0xCE5B, 0xCD51,
    0xCC4A, 0xCB48, 0xCA4A, 0xC94F, 0xC858, 0xC764, 0xC674, 0xC587, 0xC49D, 0xC3B7, 0xC2D4, 0xC1F4,
    0xC116, 0xC03C, 0xBF65, 0xBE90, 0xBDBE, 0xBCEF, 0xBC23, 0xBB59, 0xBA91, 0xB9CC, 0xB90A, 0xB84A,
    0xB78C, 0xB6D0, 0xB617, 0xB560, 0xB451, 0xB2F0, 0xB196, 0xB044, 0xAEF9, 0xADB6, 0xAC79, 0xAB43,
    0xAA14, 0xA8EB, 0xA7C8, 0xA6AA, 0xA592, 0xA480, 0xA373, 0xA26B, 0xA168, 0xA06A, 0x9F70, 0x9E7B,
    0x9D8A, 0x9C9D, 0x9BB5, 0x9AD1, 0x99F0, 0x9913, 0x983A, 0x9765, 0x9693, 0x95C4, 0x94F8, 0x9430,
    0x936B, 0x92A9, 0x91EA, 0x912E, 0x9075, 0x8FBE, 0x8F0A, 0x8E59, 0x8DAA, 0x8CFE, 0x8C54, 0x8BAC,
    0x8B07, 0x8A64, 0x89C4, 0x8925, 0x8889, 0x87EE, 0x8756, 0x86C0, 0x862B, 0x8599, 0x8508, 0x8479,
    0x83EC, 0x8361, 0x82D8, 0x8250, 0x81C9, 0x8145, 0x80C2, 0x8040,
];

/// Calculates `(r, s)` such as `r ~= 1 / sqrt(m)` and `s ~= sqrt(m)`, where:
/// * `m` is in [1, 4), as 2.30 fixed point
/// * `r` is 0.32 fixed point, with `|r * sqrt(m) - 1| < 2^-29`
/// * `s` is 2.30 fixed point, with `|s / sqrt(m) - 1| < 2^-27.9`
#[inline]
pub(crate) fn rsqrt_sqrt_32(m: u32) -> (u32, u32) {
    // 3 as 2.30 fixed point
    const THREE: u32 = 0xC000_0000;

    // Initial approximation from a table, indexed by whether m >= 2 and
    // the 6 highest bits of the fractional part of m (or m / 2)
    let odd = m >> 31;
    let i = (odd << 6) | ((m >> (24 + odd)) & 63);
    let r = u32::from(RSQRT_TBL[i as usize]) << 16;

    // Two Goldschmidt iterations, each one roughly doubles the number of
    // correct bits:
    // s ~= m * r ~= sqrt(m)
    // d ~= s * r ~= m * r^2 ~= 1
    // u = 3 - d
    // r' = r * u / 2
    // s' = s * u / 2
    let s = mul32(m, r);
    let d = mul32(s, r);
    let u = THREE - d;
    let r = mul32(r, u) << 1;
    let s = mul32(s, u) << 1;

    let d = mul32(s, r);
    let u = THREE - d;
    let r = mul32(r, u) << 1;
    let s = mul32(s, u) << 1;

    (r, s)
}

/// Returns the high 32 bits of `a * b`
#[inline]
fn mul32(a: u32, b: u32) -> u32 {
    ((u64::from(a) * u64::from(b)) >> 32) as u32
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::Float;

    fn test<F: Float + FloatMath>() {
        use crate::sqrt;

        assert_is_nan!(sqrt(F::NAN));
        assert_is_nan!(sqrt(F::NEG_INFINITY));
        assert_is_nan!(sqrt(-F::ONE));
        assert_total_eq!(sqrt(F::INFINITY), F::INFINITY);
        assert_total_eq!(sqrt(F::ZERO), F::ZERO);
        assert_total_eq!(sqrt(-F::ZERO), -F::ZERO);
    }

    #[test]
    fn test_f32() {
        test::<f32>();
    }

    #[test]
    fn test_f64() {
        test::<f64>();
    }
}
