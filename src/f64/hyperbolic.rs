//! Hyperbolic functions for `f64`.
//!
//! They use the argument reduction and the evaluation of `exp` (see
//! `super::exp`) with `|x|` (the sign is restored at the end):
//!
//! * `sinh(x) = (exp(x) - exp(-x)) / 2`
//! * `cosh(x) = (exp(x) + exp(-x)) / 2`
//! * `tanh(x) = -q / (2 + q)`, with `q = exp(-2 * x) - 1`
//!
//! `exp(x) / 2` and `exp(-x) / 2` are evaluated as sums of two `f64` from the
//! reduced argument of `x` and its negation (`x = m * ln(2) / N + r`), and
//! added rounding only once. When `x >= 32`, `exp(-x)` is negligible.
//!
//! `sinh(x)` is much smaller than `exp(x)` when `x` is small, so the more
//! accurate evaluation of `exp_m1` is used when `x < ~ln(2)` (`k = 0`), and,
//! when `x < 2^-9` (so `m = 0`), `sinh(x)` and `cosh(x)` are calculated
//! directly as the odd and even parts of the polynomial of `exp_m1`.

use super::exp::{
    EXP2_TBL_BITS, Reduced, exp_m1_parts, exp_m1_poly_parts, exp_poly, exp_poly_parts,
};
use super::f64x2::F64x2;
use crate::traits::Float as _;

impl crate::generic::Hyperbolic for f64 {
    #[inline]
    fn sinh_finite(x: Self) -> Self {
        let x_exp = x.exponent();
        if x_exp < -26 {
            // sinh(x) = x + x^3 / 6 + ..., where |x| < 2^-26, so
            // x^3 / 6 + ... < 2^-54 * |x|, which is less than half the spacing
            // of `f64` around x, and the result rounds to x
            x
        } else if x_exp >= 10 {
            f64::INFINITY.copysign(x)
        } else {
            sinh_cosh_abs(x.abs()).0.copysign(x)
        }
    }

    #[inline]
    fn cosh_finite(x: Self) -> Self {
        let x_exp = x.exponent();
        if x_exp < -26 {
            // cosh(x) = 1 + x^2 / 2 + ..., where |x| < 2^-26, so
            // x^2 / 2 + ... < 2^-53, and the result rounds to 1
            1.0
        } else if x_exp >= 10 {
            f64::INFINITY
        } else {
            cosh_abs(x.abs())
        }
    }

    #[inline]
    fn sinh_cosh_finite(x: Self) -> (Self, Self) {
        let x_exp = x.exponent();
        if x_exp < -26 {
            // see `sinh_finite` and `cosh_finite`
            (x, 1.0)
        } else if x_exp >= 10 {
            (f64::INFINITY.copysign(x), f64::INFINITY)
        } else {
            let (sinh, cosh) = sinh_cosh_abs(x.abs());
            (sinh.copysign(x), cosh)
        }
    }

    #[inline]
    fn tanh_finite(x: Self) -> Self {
        let x_exp = x.exponent();
        if x_exp < -27 {
            // tanh(x) = x - x^3 / 3 + ..., where |x| < 2^-27, so
            // x^3 / 3 - ... < 2^-55 * |x|, which is less than half the spacing
            // of `f64` around x, and the result rounds to x
            x
        } else if x_exp >= 5 {
            // |x| >= 32, 1 - tanh(|x|) < 2^-91
            1.0f64.copysign(x)
        } else {
            tanh_abs(x.abs()).copysign(x)
        }
    }
}

/// Returns `(sinh(x), cosh(x))`, for `2^-26 <= x < 1024`.
///
/// It is always inlined, so the unused result is not calculated.
#[inline(always)]
fn sinh_cosh_abs(x: f64) -> (f64, f64) {
    if x.exponent() < -9 {
        // x < 2^-9, see `small_sinh_cosh`
        return small_sinh_cosh(x);
    }

    let red = Reduced::exp(x);
    if x < 32.0 {
        // exp(x) / 2 = ph + pl
        // exp(-x) / 2 = nh + nl
        let pos = red.mul_exp2(-1);
        let neg = red.neg().mul_exp2(-1);
        let ((ph, pl), (nh, nl)) = if (red.m >> EXP2_TBL_BITS) == 0 {
            // x < ~ln(2), sinh(x) can be much smaller than exp(x)
            // Q(±r) = qa ± r * qb
            let (qa, qb) = exp_m1_poly_parts(red.r);
            let rqb = red.r * qb;
            (pos.eval_precise(qa + rqb), neg.eval_precise(qa - rqb))
        } else {
            // x > ~ln(2), sinh(x) > exp(x) / 2.7
            // P(±r) = pa ± r * pb
            let (pa, pb) = exp_poly_parts(red.r);
            let rpb = red.r * pb;
            (pos.eval(pa + rpb), neg.eval(pa - rpb))
        };

        // exp(x) > exp(-x), the sums are exact
        let s = F64x2::fast_add11(ph, -nh);
        let c = F64x2::fast_add11(ph, nh);
        let sinh = s.hi() + (s.lo() + (pl - nl));
        let cosh = c.hi() + (c.lo() + (pl + nl));
        (sinh, cosh)
    } else {
        // exp(-x) < 2^-92 * exp(x), so sinh(x) and cosh(x) round to the same
        // value as exp(x) / 2, which can overflow
        let y = red.mul_exp2(-1).eval_to_f64(exp_poly(red.r));
        (y, y)
    }
}

/// Returns `cosh(x)`, for `2^-26 <= x < 1024`.
///
/// Unlike `sinh`, there is no cancellation, so it does not need the more
/// accurate evaluations of `sinh_cosh_abs`.
#[inline]
fn cosh_abs(x: f64) -> f64 {
    if x.exponent() < -9 {
        // x < 2^-9, see `small_sinh_cosh`
        return small_sinh_cosh(x).1;
    }

    let red = Reduced::exp(x);
    if x < 32.0 {
        // exp(x) / 2 = ph + pl
        // exp(-x) / 2 = nh + nl
        // P(±r) = pa ± r * pb
        let (pa, pb) = exp_poly_parts(red.r);
        let rpb = red.r * pb;
        let (ph, pl) = red.mul_exp2(-1).eval(pa + rpb);
        let (nh, nl) = red.neg().mul_exp2(-1).eval(pa - rpb);

        // exp(x) >= exp(-x), the sum is exact
        let c = F64x2::fast_add11(ph, nh);
        c.hi() + (c.lo() + (pl + nl))
    } else {
        // see `sinh_cosh_abs`
        red.mul_exp2(-1).eval_to_f64(exp_poly(red.r))
    }
}

/// Returns `(sinh(x), cosh(x))`, for `0 <= x < 2^-9`.
///
/// It is always inlined, so the unused result is not calculated.
#[inline(always)]
fn small_sinh_cosh(x: f64) -> (f64, f64) {
    // exp(x) = 1 + x + x^2 / 2 + x^3 * Q(x), where `x^3 * Q(x)` has an even
    // part `x^4 * qb` and an odd part `x^3 * qa` (see `exp_m1_poly_parts`)
    // sinh(x) = x + x^3 * qa
    // cosh(x) = 1 + x^2 * (1/2 + x^2 * qb)
    // where the errors of the terms after `x` and `1` (from the polynomial and
    // the roundings) are less than about 2^-68 relative to the results
    let (qa, qb) = exp_m1_poly_parts(x);
    let x2 = x * x;
    (x + (x2 * x) * qa, 1.0 + x2 * (0.5 + x2 * qb))
}

/// Returns `tanh(x)`, for `2^-27 <= x < 32`.
#[inline]
fn tanh_abs(x: f64) -> f64 {
    // tanh(x) = -q / (2 + q)
    // q = exp(-2 * x) - 1 = qh + ql
    let (qh, ql) = exp_m1_parts(-2.0 * x);

    // d = 2 + q = dh + dl, with |qh| < 1
    let d = F64x2::fast_add11(2.0, qh);
    let (dh, dl) = (d.hi(), d.lo() + ql);
    let inv_d = 1.0 / (dh + dl);

    // y0 ~= -q / d, with 26 bits, so the residual of the division can be
    // calculated exactly with a 26-bit split of `dh`
    let y0 = (-(qh + ql) * inv_d).split_hi();
    let dh_hi = dh.split_hi();
    let dh_lo = (dh - dh_hi) + dl;

    // residual = -q - y0 * d
    // `y0 * dh_hi` is exact (26 bits * 26 bits), and so is `-qh - y0 * dh_hi`
    // (by Sterbenz lemma), because `y0 * dh_hi = -qh * (1 + e)`, with `|e|`
    // of at most about 2^-12: `|ql|` is up to about 2^-13 * |qh| (in the
    // small-argument path of `exp_m1_parts`), and the truncations of `y0` and
    // `dh_hi` add less than 2^-24
    let res = ((-qh - y0 * dh_hi) - ql) - y0 * dh_lo;
    y0 + res * inv_d
}
