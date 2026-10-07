//! Exponential functions for `f32`.
//!
//! They use the table of the `f64` functions (see `crate::f64::exp`), but the
//! evaluation uses plain `f64` arithmetic. The argument is converted to
//! `z = x * log2(B) * N`, where `B` is the base and `N = 2^EXP2_TBL_BITS`,
//! which is reduced to `z = m + f`, with an integer `m` and `|f| <= ~1/2`, so
//!
//! `B^x = 2^(m / N) * 2^(f / N)`
//!
//! `S = 2^(m / N)` is taken from the table (rounded to `f64`) and
//! `2^(f / N) - 1 = exp(r) - 1`, with `r = f * ln(2) / N`, is approximated
//! with a polynomial in `f`. The coefficients are calculated from those of a
//! polynomial in `r`.

use crate::f64::{EXP2_TBL_BITS, exp2_tbl, round_i32};
use crate::traits::Float as _;

// GENERATE: consts f64 LOG2_E LOG2_10 LN_2
const LOG2_E: f64 = f64::from_bits(0x3FF71547652B82FE); // 1.4426950408889634e0
const LOG2_10: f64 = f64::from_bits(0x400A934F0979A371); // 3.321928094887362e0
const LN_2: f64 = f64::from_bits(0x3FE62E42FEFA39EF); // 6.931471805599453e-1

/// `N`
const TBL_N: f64 = (1 << EXP2_TBL_BITS) as f64;

/// `ln(2) / N`
const LN_2_N: f64 = LN_2 / TBL_N;

impl crate::generic::Exp for f32 {
    #[inline]
    fn exp_finite(x: Self) -> Self {
        if x.exponent() >= 7 {
            if x.is_sign_negative() {
                return 0.0;
            } else {
                return f32::INFINITY;
            }
        }

        exp_f64(f64::from(x)) as f32
    }

    #[inline]
    fn exp_m1_finite(x: Self) -> Self {
        if x.exponent() >= 7 {
            if x.is_sign_negative() {
                return -1.0;
            } else {
                return f32::INFINITY;
            }
        }

        // exp(x) - 1 = (S - 1) + S * (2^(f / N) - 1)
        //
        // The relative error before the final rounding is about 2^-43.5:
        // * m = 0: S = 1, so the result is the polynomial.
        // * m != 0: |exp(x) - 1| > 0.0027 * S, and the polynomial has an
        //   absolute error of 2^-52.1 relative to `S`. `S - 1` has an error of
        //   at most 2^-53 * |S - 1| (it is exact when 1/2 <= S <= 2^53), and
        //   the error of `x * log2(e) * N` is about 2^-53 * |x| relative to
        //   `S`, which is also small compared to |exp(x) - 1| / S.
        let (s, f) = reduce(f64::from(x) * (LOG2_E * TBL_N));
        ((s - 1.0) + (s * f) * poly_m1(f)) as f32
    }

    #[inline]
    fn exp2_finite(x: Self) -> Self {
        if x.exponent() >= 8 {
            if x.is_sign_negative() {
                return 0.0;
            } else {
                return f32::INFINITY;
            }
        }

        // Like `exp`, but `x * N` is exact
        let (s, f) = reduce(f64::from(x) * TBL_N);
        (s + (s * f) * poly(f)) as f32
    }

    #[inline]
    fn exp10_finite(x: Self) -> Self {
        if x.exponent() >= 6 {
            if x.is_sign_negative() {
                return 0.0;
            } else {
                return f32::INFINITY;
            }
        }

        // Like `exp`, `x * log2(10) * N` (|x| < 64) adds an error of 2^-46
        let (s, f) = reduce(f64::from(x) * (LOG2_10 * TBL_N));
        (s + (s * f) * poly(f)) as f32
    }
}

/// Returns `exp(x)`, for `|x| < 128`, with a relative error of about 2^-40.5:
/// * Polynomial: 2^-40.7 (relative to `exp(r) ~= 1`).
/// * `x * log2(e) * N`: 2^-46.
/// * Table value and evaluation: ~2^-51.
#[inline]
pub(super) fn exp_f64(x: f64) -> f64 {
    let (s, f) = reduce(x * (LOG2_E * TBL_N));
    s + (s * f) * poly(f)
}

/// Splits `z = m + f`, where `m` is an integer and `|f| <= ~1/2`, returning
/// `(2^(m / N), f)`. `2^(m / N)` is rounded to `f64`.
///
/// `|z|` must be less than 2^15.
#[inline]
pub(super) fn reduce(z: f64) -> (f64, f64) {
    let (mf, m) = round_i32(z);
    let (s, _) = exp2_tbl(m);
    // exact
    (s, z - mf)
}

/// Returns `(2^(f / N) - 1) / f`, with a relative error of 2^-32.2 (in
/// `2^(f / N) - 1`) for `|f| <= ~1/2`.
#[inline]
fn poly(f: f64) -> f64 {
    // exp(r) - 1 ~= r + K2 * r^2 + K3 * r^3 for |r| <= 0.002711
    // GENERATE: exp_m1_poly f64 2 -0.002711 0.002711
    const K2: f64 = f64::from_bits(0x3FE000007B4B2257); // 5.000002296520262e-1
    const K3: f64 = f64::from_bits(0x3FC55555D8DBB03A); // 1.6666672791268572e-1

    // The same polynomial, with `r = f * ln(2) / N`
    const P1: f64 = LN_2_N;
    const P2: f64 = K2 * (LN_2_N * LN_2_N);
    const P3: f64 = K3 * (LN_2_N * LN_2_N * LN_2_N);

    P1 + f * (P2 + f * P3)
}

/// `[P1, P2, P3, P4]` such that
/// `2^(f / N) - 1 ~= f * (P1 + P2 * f + P3 * f^2 + P4 * f^3)`, with a relative
/// error of 2^-43.6 for `|f| <= ~1/2`.
const EXP2_M1_POLY: [f64; 4] = {
    // exp(r) - 1 ~= r + K2 * r^2 + K3 * r^3 + K4 * r^4 for |r| <= 0.002711 (the
    // range includes a margin for the double rounding of x87 in `round_i32`)
    // GENERATE: exp_m1_poly f64 3 -0.002711 0.002711
    const K2: f64 = f64::from_bits(0x3FDFFFFFFFFFFDD0); // 4.999999999999689e-1
    const K3: f64 = f64::from_bits(0x3FC55555C24AC4A7); // 1.6666671740452907e-1
    const K4: f64 = f64::from_bits(0x3FA55555D1560E5C); // 4.1666681102495245e-2

    // The same polynomial, with `r = f * ln(2) / N`
    [
        LN_2_N,
        K2 * (LN_2_N * LN_2_N),
        K3 * (LN_2_N * LN_2_N * LN_2_N),
        K4 * (LN_2_N * LN_2_N * LN_2_N * LN_2_N),
    ]
};

/// Returns `(2^(f / N) - 1) / f` (see `EXP2_M1_POLY`).
#[inline]
fn poly_m1(f: f64) -> f64 {
    let [p1, p2, p3, p4] = EXP2_M1_POLY;
    p1 + f * (p2 + f * (p3 + f * p4))
}

/// Returns `(even, odd)`, the even and odd parts of `2^(f / N) - 1` (see
/// `EXP2_M1_POLY`), so `2^(±f / N) - 1 ~= even ± odd`.
#[inline]
pub(super) fn exp2_m1_parts(f: f64) -> (f64, f64) {
    let [p1, p2, p3, p4] = EXP2_M1_POLY;
    let f2 = f * f;
    (f2 * (p2 + f2 * p4), f * (p1 + f2 * p3))
}
