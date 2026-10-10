//! Power functions for `f64`.
//!
//! `|x|^y = exp(y * ln(|x|))`, where `ln(|x|)` is calculated with extra
//! precision (see `super::log::ln_accurate_parts`), because its relative error
//! is amplified by `|y * ln(|x|)|`, which can be up to about 745 when the
//! result is finite and non-zero. `y * ln(|x|)` is calculated as a sum of two
//! `f64`, and its exponential with the evaluation of `exp` (see `super::exp`).
//!
//! The relative error before the final rounding is about 2^-64.7:
//! * Logarithm: 2^-75 (see `super::log::ln_accurate_parts`), amplified by up
//!   to 2^9.6.
//! * Product: 2^-76 (see `exp_y_ln`), also amplified by up to 2^9.6.
//! * Exponential: 2^-68 (see `super::exp::Reduced::eval`).
//!
//! Integer exponents with `|y| <= 64` (in both `pow` and `powi`) use binary
//! exponentiation with extra precision instead (see `powi_small`), when the
//! result is far from overflowing or underflowing.

use super::exp::{Reduced, exp_poly};
use super::f64x2::F64x2;
use super::log::ln_accurate_parts;
use crate::generic::is_int;
use crate::traits::Float as _;

impl crate::generic::Pow for f64 {
    #[inline]
    fn pow_finite(x: Self, xedelta: Self::Exp, y: Self, sign: bool) -> Self {
        // Small integer exponents are evaluated like in `powi`, which is
        // faster and gives the same results
        let absz = if y.abs() <= f64::from(POWI_SMALL_MAX) && is_int(y) {
            abs_powi(x, xedelta, y as i32)
        } else {
            exp_y_ln(x.abs(), xedelta, y)
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
fn abs_powi(x: f64, xedelta: i16, y: i32) -> f64 {
    let n = y.unsigned_abs();
    // |x| = 2^e * m, with 1 <= m < 2, so 2^(n * e) <= |x^n| < 2^(n * (e + 1))
    let e = i32::from(x.exponent());
    if xedelta == 0
        && n <= POWI_SMALL_MAX
        && (n as i32) * e >= -POWI_SMALL_RANGE
        && (n as i32) * (e + 1) <= POWI_SMALL_RANGE
    {
        powi_small(x.abs(), y)
    } else {
        exp_y_ln(x.abs(), xedelta, f64::from(y))
    }
}

/// Maximum `|y|` for `powi_small`.
const POWI_SMALL_MAX: u32 = 64;

/// `x^y` and its reciprocal must be in `[2^-RANGE, 2^RANGE]` for `powi_small`,
/// so no intermediate value overflows or is subnormal.
const POWI_SMALL_RANGE: i32 = 960;

/// Returns `x^y`, for `1 <= |y| <= POWI_SMALL_MAX`, with binary
/// exponentiation, where `x` is positive and normal, and `x^|y|` is in
/// `[2^-POWI_SMALL_RANGE, 2^POWI_SMALL_RANGE]`.
///
/// The powers `x^m` (with `m <= |y|`, so they are between 1 and `x^|y|`)
/// are calculated as sums of two `f64`, `h + l`, where `h` has 26 bits and
/// `|l| < 2^-24.9 * h`, so the leading products `h * h` and `h * x1` (where
/// `x = x1 + x2` and `x1` has 26 bits) are exact. Each squaring or product
/// has a relative error less than about 2^-75.6, so `x^|y|` has a relative
/// error less than `(|y| - 1) * 2^-75.6`, which is 2^-69.6 at most.
#[inline]
fn powi_small(x: f64, y: i32) -> f64 {
    let n = y.unsigned_abs();
    let x1 = x.split_hi();
    let x2 = x - x1;

    // From the most significant bit of `n`
    let (mut h, mut l) = (x1, x2);
    let mut bit = (1u32 << (31 - n.leading_zeros())) >> 1;
    while bit != 0 {
        // (h + l)^2 = h^2 + (2 * h + l) * l
        (h, l) = renorm(h * h, (2.0 * h + l) * l);
        if (n & bit) != 0 {
            // (h + l) * x = h * x1 + (h * x2 + l * x)
            (h, l) = renorm(h * x1, h * x2 + l * x);
        }
        bit >>= 1;
    }

    if y < 0 {
        // 1 / (h + l) = q / (1 - r) ~= q + q * r, where `q ~= 1 / (h + l)`
        // and r = 1 - q * (h + l), with |r| < 2^-51
        //
        // q = q1 + q2, where `q1` has 26 bits, so `q1 * h` and `q2 * h` are
        // exact, and so is `1 - q1 * h` (by Sterbenz lemma). The remaining
        // terms are less than 2^-24 in magnitude, so the error of `r` is less
        // than about 2^-76.
        let q = (1.0 / (h + l)).purify();
        let q1 = q.split_hi();
        let q2 = q - q1;
        let r = ((1.0 - q1 * h) - q2 * h) - q * l;
        q + q * r
    } else {
        h + l
    }
}

/// Renormalizes `p + m` to `h + l`, where `h` has 26 bits and
/// `|l| < 2^-24.9 * |h|`, for `|m| <= 2^-22 * |p|`, with a relative error less
/// than 2^-77.9.
#[inline]
fn renorm(p: f64, m: f64) -> (f64, f64) {
    // `h` is `p + m` rounded and truncated to 26 bits, so it is a multiple of
    // the ULP of `p`, and `p - h` is exact (its magnitude is less than
    // 2^-20 * |p|).
    let h = (p + m).purify().split_hi();
    (h, (p - h) + m)
}

/// Returns `exp(y * ln(x * 2^xedelta))`, where `x` is positive and normal.
#[inline]
fn exp_y_ln(x: f64, xedelta: i16, y: f64) -> f64 {
    // ln(x) = l_hi + l_lo
    let (l_hi, l_lo) = ln_accurate_parts(x, xedelta);

    // y * ln(x) = p_hi + p_lo
    //
    // ln(x) = a + b and y = y1 + y2, where `a` (`l_hi` truncated) and `y1`
    // have 26 bits, so `y1 * a` and `y2 * a` are exact, and
    // y * ln(x) = y1 * a + (y2 * a + y * b)
    // The terms in parentheses are less than 2^-24 * |y * ln(x)|, so their
    // rounding errors (and the one of `b`) are less than 2^-76 * |y * ln(x)|.
    let a = l_hi.split_hi();
    let b = (l_hi - a) + l_lo;
    let y1 = y.split_hi();
    let y2 = y - y1;
    let p = F64x2::fast_add11(y1 * a, y2 * a + y * b);

    if p.hi().exponent() >= 10 {
        // |y * ln(x)| >= 1024 (or it overflows), so the result overflows or
        // underflows
        if p.hi().is_sign_negative() {
            0.0
        } else {
            f64::INFINITY
        }
    } else {
        let red = Reduced::exp_sum(p.hi(), p.lo());
        red.eval_to_f64(exp_poly(red.r))
    }
}
