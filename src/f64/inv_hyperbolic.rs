//! Inverse hyperbolic functions for `f64`.
//!
//! They use the evaluation of `ln` (see `super::log`) with `|x|` (the sign is
//! restored at the end):
//!
//! * `asinh(x) = ln(x + sqrt(x^2 + 1))`
//! * `acosh(x) = ln(x + sqrt(x^2 - 1))`
//! * `atanh(x) = ln((1 + x) / (1 - x)) / 2`
//!
//! The arguments of `ln` are calculated as sums of two `f64`, with relative
//! errors that are small compared to the one of `ln` (see `asinh_parts`,
//! `acosh_parts` and `ln_atanh_parts`).
//!
//! When `x >= 16`, `asinh` and `acosh` avoid the square root:
//!
//! * `asinh(x) = ln(2 * x) + ln((1 + sqrt(1 + 1 / x^2)) / 2)`
//! * `acosh(x) = ln(2 * x) + ln((1 + sqrt(1 - 1 / x^2)) / 2)`
//!
//! where the second terms are approximated with polynomials in `1 / x^2`.
//!
//! When `|x|` is small (less than 2^-5 for `asinh` and 2^-6 for `atanh`),
//! `asinh` and `atanh` are approximated with polynomials.

use super::f64x2::F64x2;
use super::log::{ln_parts, ln_sum_parts};
use super::{sqrt_parts, square_parts};
use crate::traits::Float as _;

impl crate::generic::InvHyperbolic for f64 {
    #[inline]
    fn asinh_finite(x: Self) -> Self {
        let x_exp = x.exponent();
        if x_exp < -26 {
            // asinh(x) = x - x^3 / 6 + ..., where |x| < 2^-26, so
            // x^3 / 6 - ... < 2^-54 * |x|, which is less than half the spacing
            // of `f64` around x, and the result rounds to x
            x
        } else if x_exp < -5 {
            // |x| < 2^-5
            asinh_small(x)
        } else {
            // asinh(x) = asinh(|x|) * sgn(x)
            let (hi, lo) = asinh_parts(x.abs());
            (hi + lo).copysign(x)
        }
    }

    #[inline]
    fn acosh_finite(x: Self) -> Self {
        let (hi, lo) = acosh_parts(x);
        hi + lo
    }

    #[inline]
    fn atanh_finite(x: Self) -> Self {
        let x_exp = x.exponent();
        if x_exp < -27 {
            // atanh(x) = x + x^3 / 3 + ..., where |x| < 2^-27, so
            // x^3 / 3 + ... < 2^-54 * |x|, which is less than half the spacing
            // of `f64` around x, and the result rounds to x
            x
        } else if x_exp < -6 {
            // |x| < 2^-6
            atanh_small(x)
        } else {
            // atanh(x) = ln((1 + |x|) / (1 - |x|)) / 2 * sgn(x)
            let (hi, lo) = ln_atanh_parts(x.abs());
            (0.5 * (hi + lo)).copysign(x)
        }
    }
}

/// Returns `asinh(x)`, for `2^-26 <= |x| < 2^-5`.
#[inline]
fn asinh_small(x: f64) -> f64 {
    // asinh(x) ~= x + x^3 * (K3 + K5 * x^2 + ... + K11 * x^8)
    // GENERATE: asinh_poly f64 5 -0.0313 0.0313
    const K3: f64 = f64::from_bits(0xBFC5555555555555); // -1.6666666666666666e-1
    const K5: f64 = f64::from_bits(0x3FB33333333331A1); // 7.499999999999442e-2
    const K7: f64 = f64::from_bits(0xBFA6DB6DB6A1CF5F); // -4.464285711665527e-2
    const K9: f64 = f64::from_bits(0x3F9F1C6E2BE455B5); // 3.0381890706015564e-2
    const K11: f64 = f64::from_bits(0xBF96DB9BE588322A); // -2.2322116741678645e-2

    // -1/6 = K3_HI + K3_LO, where `K3_HI` (`K3` truncated) has 11 bits and
    // `K3_LO` has a relative error of 2^-53 (2^-65 relative to 1/6). The
    // rounding of `K3` alone is 2^-66.6 relative to the result.
    const K3_HI: f64 = f64::from_bits(K3.to_bits() & (u64::MAX << 42));
    const K3_LO: f64 = (-1.0 - 6.0 * K3_HI) / 6.0;

    // The cubic term is up to 2^-12.5 * |x|, so its rounding errors would be
    // up to about 2^-64 relative to the result. Its leading part is added
    // exactly to `x` instead: x = a + b, where `a` has 14 bits, so
    // `K3_HI * a^3` is exact, and
    // -x^3 / 6 = K3_HI * a^3 + K3_HI * b * (x^2 + x * a + a^2) + K3_LO * x^3
    // The remaining terms are less than 2^-22 * |x|, so their rounding errors
    // are less than about 2^-73 relative to the result.
    let a = f64::from_bits(x.to_bits() & (u64::MAX << 39));
    let b = x - a;
    let a2 = a * a;
    let s = F64x2::fast_add11(x, K3_HI * (a2 * a));

    let x2 = x * x;
    let x4 = x2 * x2;
    let c3 = K3_HI * (b * ((x2 + x * a) + a2)) + K3_LO * (x2 * x);
    let q = (K5 + x2 * K7) + x4 * (K9 + x2 * K11);
    s.hi() + ((s.lo() + c3) + (x4 * x) * q)
}

/// Returns `atanh(x)`, for `2^-27 <= |x| < 2^-6`.
#[inline]
fn atanh_small(x: f64) -> f64 {
    // atanh(x) ~= x + x^3 * (K3 + K5 * x^2 + K7 * x^4 + K9 * x^6)
    // GENERATE: atanh_poly f64 4 -0.0157 0.0157
    const K3: f64 = f64::from_bits(0x3FD5555555555555); // 3.333333333333333e-1
    const K5: f64 = f64::from_bits(0x3FC99999999A1B97); // 2.0000000000092363e-1
    const K7: f64 = f64::from_bits(0x3FC2492479FC2D72); // 1.4285713154127283e-1
    const K9: f64 = f64::from_bits(0x3FBC755BB73FF572); // 1.1116574500916501e-1

    // The terms after `x` are less than 2^-13.5 * |x|, so their rounding
    // errors are less than about 2^-65 relative to the result.
    let x2 = x * x;
    let x4 = x2 * x2;
    let q = (K3 + x2 * K5) + x4 * (K7 + x2 * K9);
    x + x * (x2 * q)
}

/// Returns `(hi, lo)` such that `hi + lo ~= asinh(x)`, for `x >= 2^-5`.
#[inline]
fn asinh_parts(x: f64) -> (f64, f64) {
    if x.exponent() >= 4 {
        // x >= 16
        ln_2x_plus(x, 0.25, asinh_large_tail)
    } else {
        // w = x^2 + 1 = w_hi + w_lo, where the sum is exact and the error of
        // `square_parts` is negligible
        let (p, pe) = square_parts(x);
        let w = F64x2::add11(1.0, p);
        // y = sqrt(x^2 + 1) = y_hi + y_lo, with the relative error of
        // `sqrt_parts`
        let (y_hi, y_lo) = sqrt_parts(w.hi(), w.lo() + pe);
        // t = x + y = t_hi + t_lo, where y > x, so the relative error of `t`
        // is at most about the one of `y`
        let t = F64x2::fast_add11(y_hi, x);
        // asinh(x) = ln(t)
        ln_sum_parts(t.hi(), t.lo() + y_lo)
    }
}

/// Returns `(hi, lo)` such that `hi + lo ~= acosh(x)`, for `x > 1`.
#[inline]
fn acosh_parts(x: f64) -> (f64, f64) {
    if x.exponent() >= 4 {
        // x >= 16
        ln_2x_plus(x, -0.25, acosh_large_tail)
    } else {
        // w = x^2 - 1 = w_hi + w_lo, which must be exact when `x` is close
        // to one, where `w` is small: `p + pe` is exact when
        // `x < 1 + 2^-26`, and so are `p - 1` and `w_lo + pe` (when `p <= 2`).
        // Otherwise, w > 2^-25, so the errors are less than 2^-77 * w.
        let (p, pe) = square_parts(x);
        let w = F64x2::fast_add11(p, -1.0);
        // `w_hi >= 2^-51`, but `|pe|` can be up to 2^-53, so renormalize to
        // keep `w_lo` much smaller than `w_hi`, as `sqrt_parts` requires
        let w = F64x2::fast_add11(w.hi(), w.lo() + pe);
        // y = sqrt(x^2 - 1) = y_hi + y_lo, with half the relative error of
        // `w`, plus the one of `sqrt_parts`
        let (y_hi, y_lo) = sqrt_parts(w.hi(), w.lo());
        // t = x + y = t_hi + t_lo, where y < x. The error of `w` adds at
        // most about 2^-90 relative to `t` (when `x` is close to 1 + 2^-26,
        // where it is the largest relative to `w`, and y ~= 2^-12.5), which is
        // more than the error of `sqrt_parts`.
        let t = F64x2::fast_add11(x, y_hi);
        // acosh(x) = ln(t), which is greater than 2^-26 (`x >= 1 + 2^-52`),
        // as `ln_sum_parts` requires
        ln_sum_parts(t.hi(), t.lo() + y_lo)
    }
}

/// Returns `(hi, lo)` such that `hi + lo ~= 2 * atanh(x)`, for
/// `2^-6 <= x < 1`.
#[inline]
fn ln_atanh_parts(x: f64) -> (f64, f64) {
    // u = 1 + x = u_hi + u_lo
    // v = 1 - x = v_hi + v_lo
    let u = F64x2::fast_add11(1.0, x);
    let v = F64x2::fast_add11(1.0, -x);

    // t = u / v = q + r / v, where `q ~= u_hi / v_hi` and r = u - q * v
    //
    // q = q1 + q2 and v = v1 + v2, where `q1` and `v1` (`v_hi` truncated)
    // have 26 bits, so `q1 * v1` and `q2 * v1` are exact, and so is
    // `u_hi - q1 * v1` (by Sterbenz lemma). The remaining terms (including
    // `v2`) are less than 2^-23 * u_hi, so their rounding errors are less than
    // about 2^-75 * u_hi, and so relative to `t`, which is less than 2^-70
    // relative to the result (`ln(t)` is greater than 2^-5).
    let q = (u.hi() / v.hi()).purify();
    let q1 = q.split_hi();
    let q2 = q - q1;
    let v1 = v.hi().split_hi();
    let v2 = (v.hi() - v1) + v.lo();
    let r = (((u.hi() - q1 * v1) - q2 * v1) + u.lo()) - q * v2;

    // 2 * atanh(x) = ln(t)
    ln_sum_parts(q, r / v.hi())
}

/// Returns `(hi, lo)` such that `hi + lo ~= ln(2 * x) + g(1 / x^2)`, for
/// `x >= 16`, where `g(u) = k * u + tail(u)` is `asinh_large_corr` (with
/// `k = 1/4` and `asinh_large_tail`) or `acosh_large_corr` (with `k = -1/4`
/// and `acosh_large_tail`).
///
/// The absolute error added to the one of `ln_parts` is less than about
/// 2^-64, while the result is greater than 3.4.
#[inline]
fn ln_2x_plus(x: f64, k: f64, tail: impl FnOnce(f64) -> f64) -> (f64, f64) {
    let (hi, lo) = ln_parts(x, 1);
    if x.exponent() >= 32 {
        // |g(1 / x^2)| < 2^-66, which is less than 2^-70 relative to the
        // result (greater than 22), so it is negligible
        (hi, lo)
    } else {
        // `k * u` is up to 2^-10, so the rounding errors of `x^2`, of the
        // division and of the sum `lo + g(u)` would add up to about 2^-61.5.
        // Instead:
        // * x^2 = p + pe (see `square_parts`), so u = u0 + du, where
        //   `u0 = 1 / p` is rounded (its error is at most 2^-62, which is
        //   2^-64 in `k * u`) and `du ~= -pe * u0^2`.
        // * `k * u0` is exact, and it is added exactly to `hi`.
        // The rounding errors of the other terms (`lo` is up to about
        // 2^-15 * hi, and `k * du + tail(u)` is less than 2^-19) are
        // negligible.
        let (p, pe) = square_parts(x);
        let u0 = (1.0 / p).purify();
        let du = -(pe * u0) * u0;
        let s = F64x2::fast_add11(hi, k * u0);
        (s.hi(), (s.lo() + lo) + (k * du + tail(u0)))
    }
}

/// Returns `ln((1 + sqrt(1 + u)) / 2)`, for `2^-500 <= u <= 2^-8` (the lower
/// bound avoids subnormal intermediate values).
#[inline]
pub(crate) fn asinh_large_corr(u: f64) -> f64 {
    0.25 * u + asinh_large_tail(u)
}

/// Returns `ln((1 + sqrt(1 + u)) / 2) - u / 4`, for `2^-500 <= u <= 2^-8`.
#[inline]
fn asinh_large_tail(u: f64) -> f64 {
    // ln((1 + sqrt(1 + u)) / 2) ~= u / 4 + u^2 * (K2 + K3 * u + ... + K6 * u^4)
    // GENERATE: asinh_acosh_large_poly f64 asinh 5 1e-30 0.00390625
    const K2: f64 = f64::from_bits(0xBFB7FFFFFFFFFFFE); // -9.374999999999997e-2
    const K3: f64 = f64::from_bits(0x3FAAAAAAAAA9ECDB); // 5.208333333299616e-2
    const K4: f64 = f64::from_bits(0xBFA17FFFFA106D8D); // -3.417968680897863e-2
    const K5: f64 = f64::from_bits(0x3F993311ECEE9C3C); // 2.4608879171548845e-2
    const K6: f64 = f64::from_bits(0xBF9319E6AC5F41DB); // -1.8653492232091878e-2

    let u2 = u * u;
    let q = (K2 + u * K3) + u2 * ((K4 + u * K5) + u2 * K6);
    u2 * q
}

/// Returns `ln((1 + sqrt(1 - u)) / 2)`, for `2^-500 <= u <= 2^-8` (the lower
/// bound avoids subnormal intermediate values).
#[inline]
pub(crate) fn acosh_large_corr(u: f64) -> f64 {
    -0.25 * u + acosh_large_tail(u)
}

/// Returns `ln((1 + sqrt(1 - u)) / 2) + u / 4`, for `2^-500 <= u <= 2^-8`.
#[inline]
fn acosh_large_tail(u: f64) -> f64 {
    // ln((1 + sqrt(1 - u)) / 2) ~= -u / 4 + u^2 * (K2 + K3 * u + ... + K6 * u^4)
    // GENERATE: asinh_acosh_large_poly f64 acosh 5 1e-30 0.00390625
    const K2: f64 = f64::from_bits(0xBFB8000000000002); // -9.375000000000003e-2
    const K3: f64 = f64::from_bits(0xBFAAAAAAAAA9E96C); // -5.208333333299006e-2
    const K4: f64 = f64::from_bits(0xBFA180000608DB55); // -3.417968820251952e-2
    const K5: f64 = f64::from_bits(0xBF99331170D1079B); // -2.4608871947073046e-2
    const K6: f64 = f64::from_bits(0xBF93668296B98AE9); // -1.8945732545385594e-2

    let u2 = u * u;
    let q = (K2 + u * K3) + u2 * ((K4 + u * K5) + u2 * K6);
    u2 * q
}
