//! Square root for `f32`.
//!
//! The result is correctly rounded, and it is calculated with integer
//! arithmetic, like the `f64` function (see `crate::f64::sqrt`). The argument
//! is split as `x = 4^e * m`, with an integer `e` and `1 <= m < 4`, and
//! `sqrt(m)` is approximated with a relative error less than 2^-27.9 by the
//! 32-bit Goldschmidt iterations (see `crate::generic::rsqrt_sqrt_32`).
//!
//! Rounding that approximation to 24 bits gives an integer `s` that is the
//! mantissa of the correctly rounded result or differs from it by one, which
//! is decided with the remainder `m * 2^46 - s^2`, calculated exactly.

use crate::generic::rsqrt_sqrt_32;
use crate::traits::Float as _;

impl crate::generic::Sqrt for f32 {
    #[inline]
    fn sqrt_finite(x: Self) -> Self {
        let (xn, edelta) = x.normalize_arg();

        // Split x = 4^e * m such as
        // * e is an integer
        // * 1 <= m < 4
        let k = xn.exponent() + edelta;
        let odd = (k & 1) as u32;
        let e = k >> 1;
        // m as 2.30 fixed point
        let m = xn.mant() << (7 + odd);

        // s ~= sqrt(m), as 2.30 fixed point
        // |s / sqrt(m) - 1| < 2^-27.9, which is less than 1/15 of the ULP
        // of the result.
        let (_, s) = rsqrt_sqrt_32(m);

        // s ~= sqrt(m) * 2^23, rounded to an integer
        let s = (s + (1 << 6)) >> 7;

        // `s` is the correctly rounded result or differs from it by one.
        // With n = m * 2^46, the correctly rounded result satisfies
        // (s - 1/2)^2 < n < (s + 1/2)^2, which (since n and s are integers)
        // is equivalent to -s < n - s^2 <= s.
        let rem = (i64::from(m) << 16) - i64::from(s) * i64::from(s);
        let s = s + u32::from(rem > i64::from(s)) - u32::from(rem <= -i64::from(s));

        // sqrt(x) = s * 2^(e - 23), adding `s` (which has the implicit bit)
        // to the exponent minus one.
        f32::from_bits((u32::from(f32::exp_to_raw_exp(e) - 1) << 23) + s)
    }
}
