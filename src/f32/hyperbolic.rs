//! Hyperbolic functions for `f32`.
//!
//! They use the argument reduction and evaluation of `exp` (see
//! `super::exp`) with `|x|` (the sign is restored at the end). With
//! `z = |x| * log2(e) * N = m + f`:
//!
//! * `exp(±|x|) / 2 = S± * (1 + even ± odd)`, where `S± = 2^(±m / N) / 2` and
//!   `even` and `odd` are the even and odd parts of `2^(f / N) - 1`.
//! * `sinh(|x|) = d * (1 + even) + c * odd`, with `d = S+ - S-` and
//!   `c = S+ + S-`.
//! * `cosh(|x|) = c * (1 + even) + d * odd`.
//! * `tanh(|x|) = (exp(2 * |x|) - 1) / (exp(2 * |x|) + 1)`.
//!
//! When `m = 0`, `d = 0`, so `sinh` is `odd`, with no cancellation. Otherwise,
//! `|x| > ~0.0027` and the relative error of `d` (from the rounding of the
//! table values) is at most about 2^-46.

use super::exp::{exp2_m1_parts, reduce};
use crate::f64::{EXP2_TBL_BITS, exp2_tbl, round_i32};
use crate::traits::Float as _;

// GENERATE: consts f64 LOG2_E
const LOG2_E: f64 = f64::from_bits(0x3FF71547652B82FE); // 1.4426950408889634e0

/// `N`
const TBL_N: f64 = (1 << EXP2_TBL_BITS) as f64;

impl crate::generic::Hyperbolic for f32 {
    #[inline]
    fn sinh_finite(x: Self) -> Self {
        // |x| >= 2^7
        if x.raw_exp() >= Self::EXP_OFFSET + 7 {
            f32::INFINITY.copysign(x)
        } else {
            let (sinh, _) = sinh_cosh_abs(x.abs());
            (sinh as f32).copysign(x)
        }
    }

    #[inline]
    fn cosh_finite(x: Self) -> Self {
        // |x| >= 2^7
        if x.raw_exp() >= Self::EXP_OFFSET + 7 {
            f32::INFINITY
        } else {
            let (_, cosh) = sinh_cosh_abs(x.abs());
            cosh as f32
        }
    }

    #[inline]
    fn sinh_cosh_finite(x: Self) -> (Self, Self) {
        // |x| >= 2^7
        if x.raw_exp() >= Self::EXP_OFFSET + 7 {
            (f32::INFINITY.copysign(x), f32::INFINITY)
        } else {
            let (sinh, cosh) = sinh_cosh_abs(x.abs());
            ((sinh as f32).copysign(x), cosh as f32)
        }
    }

    #[inline]
    fn tanh_finite(x: Self) -> Self {
        // |x| >= 2^4, 1 - tanh(|x|) < 2^-45
        if x.raw_exp() >= Self::EXP_OFFSET + 4 {
            1.0f32.copysign(x)
        } else {
            (tanh_abs(x.abs()) as f32).copysign(x)
        }
    }
}

/// Returns `(sinh(x), cosh(x))`, for `0 <= x < 2^7`.
///
/// The relative error before the final rounding is about 2^-43.5.
#[inline(always)]
fn sinh_cosh_abs(x: f32) -> (f64, f64) {
    let z = f64::from(x) * (LOG2_E * TBL_N);
    let (mf, m) = round_i32(z);
    // exact
    let f = z - mf;

    // `exp2_tbl(m - N)` is `2^(m / N) / 2`
    let n = 1 << EXP2_TBL_BITS;
    let (sp, _) = exp2_tbl(m - n);
    let (sn, _) = exp2_tbl(-m - n);
    let d = sp - sn;
    let c = sp + sn;

    let (even, odd) = exp2_m1_parts(f);
    (d + (d * even + c * odd), c + (c * even + d * odd))
}

/// Returns `tanh(x)`, for `0 <= x < 2^4`.
///
/// The relative error before the final rounding is like the one of `exp_m1`
/// (see `super::exp`).
#[inline]
fn tanh_abs(x: f32) -> f64 {
    // exp(2 * x) = S * 2^(f / N)
    let (s, f) = reduce(f64::from(x) * (2.0 * LOG2_E * TBL_N));
    let (even, odd) = exp2_m1_parts(f);
    // S * (2^(f / N) - 1)
    let t = s * (even + odd);
    ((s - 1.0) + t) / ((s + 1.0) + t)
}
