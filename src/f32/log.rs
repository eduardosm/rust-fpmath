//! Logarithmic functions for `f32`.
//!
//! They use the tables of the `f64` functions (see `crate::f64::log`), but
//! the evaluation uses plain `f64` arithmetic. The argument is reduced to
//! `x = 2^k * (1 + z) / s`, so
//!
//! `ln(x) = (k * ln(2) - ln(s)) + ln(1 + z)`
//!
//! where `z = m * s - 1` is exact when the mantissa `m` comes from an `f32`
//! (24 bits * 24 bits), and `ln(1 + z)` is approximated with a polynomial.
//! `log2` and `log10` multiply the result by `log2(e)` and `log10(e)`.

use crate::f64::{ln_tbl, split_ln_arg};
use crate::traits::Float as _;

// GENERATE: consts f64 LOG2_E LOG10_E
const LOG2_E: f64 = f64::from_bits(0x3FF71547652B82FE); // 1.4426950408889634e0
const LOG10_E: f64 = f64::from_bits(0x3FDBCB7B1526E50E); // 4.342944819032518e-1

impl crate::generic::Log for f32 {
    #[inline]
    fn ln_finite(x: Self, edelta: i16) -> Self {
        ln_f64(f64::from(x), edelta) as f32
    }

    #[inline]
    fn ln_1p_finite(x: Self) -> Self {
        let x_exp = x.exponent();
        let x = f64::from(x);
        if x_exp < -8 {
            // |x| < 2^-8
            ln_1p_poly(x) as f32
        } else {
            // `1 + x` is exact when |x| < 2^53, otherwise its relative error
            // is at most 2^-53. `z` is not exact (`m` can have up to 53 bits),
            // but its absolute error is at most 2^-53, while the result is at
            // least about 2^-8.
            ln_f64(1.0 + x, 0) as f32
        }
    }

    #[inline]
    fn log2_finite(x: Self, edelta: i16) -> Self {
        (ln_f64(f64::from(x), edelta) * LOG2_E) as f32
    }

    #[inline]
    fn log10_finite(x: Self, edelta: i16) -> Self {
        (ln_f64(f64::from(x), edelta) * LOG10_E) as f32
    }
}

/// Returns `ln(x * 2^edelta)`, where `x` is positive and normal.
///
/// The relative error is about 2^-36.8 when the mantissa of `x` has at most
/// 29 bits (so `z` is exact):
/// * Polynomial: 2^-36.8 (relative to `ln(1 + z)`, whose magnitude is at
///   most about the one of the result).
/// * Evaluation: ~2^-52.
///
/// With wider mantissas, the rounding of `z` adds an absolute error of up to
/// 2^-53.
#[inline]
pub(super) fn ln_f64(x: f64, edelta: i16) -> f64 {
    // x * 2^edelta = 2^k * (1 + z) / s
    let (k, m, i) = split_ln_arg(x, edelta);
    let (s, t_hi, t_lo) = ln_tbl(k, i);
    let z = m * s - 1.0;
    t_hi + (ln_1p_poly(z) + t_lo)
}

/// Returns `ln(1 + z)`, for `|z| <= 2^-8`, with a relative error of about
/// 2^-36.8.
#[inline]
fn ln_1p_poly(z: f64) -> f64 {
    // ln(1 + z) ~= z + z^2 * (K2 + K3 * z + K4 * z^2)
    // GENERATE: ln_1p_poly f64 3 -0.003907 0.003907
    const K2: f64 = f64::from_bits(0xBFDFFFFFFFFB9408); // -4.9999999998391376e-1
    const K3: f64 = f64::from_bits(0x3FD5555FF0FDF601); // 3.333358624875871e-1
    const K4: f64 = f64::from_bits(0xBFD0000F1747F547); // -2.5000359796088784e-1

    let z2 = z * z;
    z + z2 * ((K2 + z * K3) + z2 * K4)
}
