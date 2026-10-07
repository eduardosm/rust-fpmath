//! Power functions for `f32`.
//!
//! `|x|^y = exp(y * ln(|x|))`, calculated with plain `f64` arithmetic, using
//! the tables of the `f64` logarithm and the `f32` exponential (see
//! `crate::f64::log` and `super::exp`).
//!
//! The relative error of `ln(|x|)` is amplified by `|y * ln(|x|)|`, which can
//! be up to about 104 when the result is finite and non-zero, so it uses the
//! polynomial of the `f64` logarithm, which is more accurate than the one of
//! the `f32` logarithm.
//!
//! The relative error before the final rounding is about 2^-40.5:
//! * Logarithm: 2^-52, amplified by up to 2^6.7.
//! * Product: 2^-53, also amplified by up to 2^6.7.
//! * Exponential: 2^-40.5.
//!
//! Integer exponents with `|y| <= 1024` (in both `pow` and `powi`) use binary
//! exponentiation in `f64` instead (see `powi_small`), when the result is far
//! from overflowing or underflowing in `f64`.

use super::exp::exp_f64;
use crate::f64::{ln_1p_q, ln_tbl, split_ln_arg};
use crate::generic::is_int;
use crate::traits::Float as _;

impl crate::generic::Pow for f32 {
    #[inline]
    fn pow_finite(x: Self, xedelta: Self::Exp, y: Self, sign: bool) -> Self {
        // Small integer exponents are evaluated like in `powi`, which is
        // faster and gives the same results
        let absz = if y.abs() <= POWI_SMALL_MAX as f32 && is_int(y) {
            abs_powi(x, xedelta, y as i32)
        } else {
            exp_y_ln(x.abs(), xedelta, f64::from(y))
        };
        if sign { -absz } else { absz }
    }

    #[inline]
    fn powi_finite(x: Self, xedelta: Self::Exp, y: i32) -> Self {
        let absz = abs_powi(x, xedelta, y);
        if x.is_sign_negative() && (y & 1) != 0 {
            -absz
        } else {
            absz
        }
    }
}

/// Returns `|x * 2^xedelta|^y`, where `x` is normal.
#[inline]
fn abs_powi(x: f32, xedelta: i16, y: i32) -> f32 {
    let n = y.unsigned_abs();
    // |x| = 2^e * m, with 1 <= m < 2, so 2^(n * e) <= |x^n| < 2^(n * (e + 1))
    let e = i32::from(x.exponent() + xedelta);
    if n <= POWI_SMALL_MAX
        && (n as i32) * e >= -POWI_SMALL_RANGE
        && (n as i32) * (e + 1) <= POWI_SMALL_RANGE
    {
        powi_small(f64::from(x.abs()) * f64::exp2i_fast(xedelta), y) as f32
    } else {
        exp_y_ln(x.abs(), xedelta, f64::from(y))
    }
}

/// Maximum `|y|` for `powi_small`.
const POWI_SMALL_MAX: u32 = 1024;

/// `x^y` and its reciprocal must be in `[2^-RANGE, 2^RANGE]` for `powi_small`,
/// so no intermediate value overflows or is subnormal.
const POWI_SMALL_RANGE: i32 = 1000;

/// Returns `x^y`, for `1 <= |y| <= POWI_SMALL_MAX`, with binary
/// exponentiation in plain `f64` arithmetic, where `x` is positive and normal,
/// and `x^|y|` is in `[2^-POWI_SMALL_RANGE, 2^POWI_SMALL_RANGE]`.
///
/// The intermediate values are powers of `x` (or of its reciprocal) with
/// exponents up to `|y|`, so they are between 1 and `x^y`.
///
/// The powers `x^(2^j)` have relative errors of about `2^j * 2^-53` (each
/// squaring doubles the error and adds a rounding), and so does the reciprocal
/// of `x` (when `y < 0`), so the relative error is less than about
/// `2 * |y| * 2^-53`, which is 2^-42 at most.
#[inline]
fn powi_small(x: f64, y: i32) -> f64 {
    let mut b = if y < 0 { 1.0 / x } else { x };
    let mut n = y.unsigned_abs();
    let mut r = 1.0;
    loop {
        if (n & 1) != 0 {
            r *= b;
        }
        n >>= 1;
        if n == 0 {
            return r;
        }
        b *= b;
    }
}

/// Returns `exp(y * ln(x * 2^xedelta))`, where `x` is positive and normal.
#[inline]
fn exp_y_ln(x: f32, xedelta: i16, y: f64) -> f32 {
    let p = y * ln(x, xedelta);
    if p.exponent() >= 7 {
        // |y * ln(x)| >= 128, so the result overflows or underflows
        if p.is_sign_negative() {
            0.0
        } else {
            f32::INFINITY
        }
    } else {
        exp_f64(p) as f32
    }
}

/// Returns `ln(x * 2^edelta)`, where `x` is positive and normal, with a
/// relative error of about 2^-52.
#[inline]
fn ln(x: f32, edelta: i16) -> f64 {
    // x * 2^edelta = 2^k * (1 + z) / s, where `z` is exact (24 bits * 24 bits)
    let (k, m, i) = split_ln_arg(f64::from(x), edelta);
    let (s, t_hi, t_lo) = ln_tbl(k, i);
    let z = m * s - 1.0;
    // ln(1 + z) ~= z - z^2 / 2 + z^3 * Q(z), where the terms after `z` are
    // less than 2^-9 * |z|
    t_hi + ((z + (z * z) * (z * ln_1p_q(z) - 0.5)) + t_lo)
}
