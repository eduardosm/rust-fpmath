//! Square root for `f64`.
//!
//! The result is correctly rounded, and it is calculated with integer
//! arithmetic. The argument is split as `x = 4^e * m`, with an integer `e`
//! and `1 <= m < 4` (as a fixed-point number), so
//!
//! `sqrt(x) = 2^e * sqrt(m)`
//!
//! `sqrt(m)` is approximated by `sqrt_approx`: an approximation of
//! `1 / sqrt(m)` from a table is refined with two Goldschmidt iterations in
//! 32-bit arithmetic (see `crate::generic::rsqrt_sqrt_32`) and one more in
//! 64-bit arithmetic, which gives `sqrt(m)`.
//!
//! Rounding that approximation to 53 bits gives an integer `s` that is the
//! mantissa of the correctly rounded result or differs from it by one, which
//! is decided with the remainder `m * 2^104 - s^2`, calculated exactly.

use crate::generic::rsqrt_sqrt_32;
use crate::traits::Float as _;

impl crate::generic::Sqrt for f64 {
    #[inline]
    fn sqrt_finite(x: Self) -> Self {
        let (xn, edelta) = x.normalize_arg();

        // Split x = 4^e * m such as
        // * e is an integer
        // * 1 <= m < 4
        let k = xn.exponent() + edelta;
        let odd = (k & 1) as u32;
        let e = k >> 1;
        // m as 2.62 fixed point
        let m = xn.mant() << (10 + odd);

        // s ~= sqrt(m) * 2^52, rounded to an integer, with an error less
        // than 1/8 before rounding
        let s = (sqrt_approx(m) + (1 << 8)) >> 9;

        // `s` is the correctly rounded result or differs from it by one.
        // With n = m * 2^104, the correctly rounded result satisfies
        // (s - 1/2)^2 < n < (s + 1/2)^2, which (since n and s are integers)
        // is equivalent to -s < n - s^2 <= s.
        // `n - s^2` is small, so it can be calculated with wrapping
        // 64-bit arithmetic.
        let rem = (m << 42).wrapping_sub(s.wrapping_mul(s)) as i64;
        let s = s + u64::from(rem > s as i64) - u64::from(rem <= -(s as i64));

        // sqrt(x) = s * 2^(e - 52), adding `s` (which has the implicit bit)
        // to the exponent minus one.
        f64::from_bits((u64::from(f64::exp_to_raw_exp(e) - 1) << 52) + s)
    }
}

/// Returns `s ~= sqrt(m)`, where:
/// * `m` is in [1, 4), as 2.62 fixed point
/// * `s` is 3.61 fixed point, with `|s / sqrt(m) - 1| < 2^-56`
#[inline]
pub(super) fn sqrt_approx(m: u64) -> u64 {
    // r ~= 1 / sqrt(m), as 0.64 fixed point
    // |r * sqrt(m) - 1| < 2^-28.6 (2^-29 from `rsqrt_sqrt_32`, plus 2^-31
    // from the truncation of m to 32 bits)
    let (r, _) = rsqrt_sqrt_32((m >> 32) as u32);
    let r = u64::from(r) << 32;

    // One more Goldschmidt iteration with 64-bit arithmetic, after which
    // the relative error is about 1.5 times the square of the one of `r`.
    let s = mul64(m, r); // 2.62
    let d = mul64(s, r); // 2.62
    let u = (3 << 62) - d; // 2.62
    mul64(s, u) // 3.61
}

/// Returns the high 64 bits of `a * b`
#[inline]
fn mul64(a: u64, b: u64) -> u64 {
    ((u128::from(a) * u128::from(b)) >> 64) as u64
}
