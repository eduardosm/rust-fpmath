//! Hypotenuse for `f32`.
//!
//! `hypot(x, y) = sqrt(x^2 + y^2)` is correctly rounded. `x^2 + y^2` is
//! calculated in `f64` (where the squares are exact), and its square root is
//! approximated with a relative error less than 2^-34 (see `sqrt_approx`).
//! Rounding that approximation to `f32` gives the correctly rounded result,
//! unless it is very close to the midpoint between two consecutive `f32`
//! values. In that case, or when both arguments are subnormal or zero, the
//! result is calculated with integer arithmetic instead (see `hypot_exact`),
//! like the `f64` function (see `crate::f64::hypot`).

use crate::generic::rsqrt_sqrt_32;
use crate::traits::Float as _;

impl crate::generic::Hypot for f32 {
    #[inline]
    fn hypot_finite(x: Self, y: Self) -> Self {
        if (x.abs().to_bits() | y.abs().to_bits()) < f32::MIN_POSITIVE.to_bits() {
            // Both are subnormal or zero, the result might be subnormal
            return hypot_exact(x, y);
        }

        // h ~= hypot(x, y) = sqrt(x^2 + y^2)
        // x^2 and y^2 are exact, and their sum and the square root approximation
        // add a relative error less than 2^-34 (see `sqrt_approx`)
        let xd = f64::from(x);
        let yd = f64::from(y);
        let h = sqrt_approx(xd * xd + yd * yd);

        // Rounding `h` to f32 gives the correctly rounded result unless the
        // exact result and `h` are on different sides of the midpoint between
        // two f32 numbers, which can only happen when `h` is close to one.
        // Midpoints are at 2^28 in the 29 lowest bits of `h`, and a relative
        // error of 2^-32 (with some margin) is less than 2^21 units of those bits.
        const TOL: u64 = 1 << 21;
        if ((h.to_bits() & ((1 << 29) - 1)).wrapping_sub((1 << 28) - TOL)) <= 2 * TOL {
            return hypot_exact(x, y);
        }

        h as f32
    }
}

/// Calculates `sqrt(x)` for a positive and normal `x`, with a relative
/// error less than 2^-34.
#[inline]
fn sqrt_approx(x: f64) -> f64 {
    // r ~= 1 / sqrt(x), from bit fiddling, with a relative error less than 0.035
    let hx = 0.5 * x;
    let mut r = f64::from_bits(0x5FE6_EB50_C7B5_37A9_u64.wrapping_sub(x.to_bits() >> 1));
    // Three Newton iterations, the relative error is roughly squared in each
    // one: 0.035 -> 2^-9.2 -> 2^-17.7 -> 2^-34.8
    // r_next = r * (3 - x * r^2) / 2
    r *= 1.5 - hx * r * r;
    r *= 1.5 - hx * r * r;
    r *= 1.5 - hx * r * r;
    r * x
}

/// Calculates `hypot(x, y)` with exact integer arithmetic.
#[cold]
fn hypot_exact(x: f32, y: f32) -> f32 {
    // a = max(|x|, |y|) and b = min(|x|, |y|), as raw bits
    let xraw = x.abs().to_bits();
    let yraw = y.abs().to_bits();
    let (a, b) = if xraw >= yraw {
        (xraw, yraw)
    } else {
        (yraw, xraw)
    };

    if b == 0 {
        // hypot(a, 0) = a
        return f32::from_bits(a);
    }

    if a < f32::MIN_POSITIVE.to_bits() {
        // Both are subnormal, so hypot(a, b) = sqrt(A^2 + B^2) * 2^-149,
        // where A and B are the raw bits of `a` and `b`. The result is
        // less than 2^-125, so its ULP is 2^-149 and its raw bits are
        // sqrt(A^2 + B^2) rounded to an integer.
        return f32::from_bits(isqrt_round(
            u64::from(a) * u64::from(a) + u64::from(b) * u64::from(b),
        ));
    }

    // a = ma * 2^(ea - 23), with 2^23 <= ma < 2^24
    // b = mb * 2^(eb - 23), with 2^23 <= mb < 2^24
    let a = f32::from_bits(a);
    let ea = a.exponent();
    let ma = a.mant();
    let (b, edelta) = f32::from_bits(b).normalize_arg();
    let eb = b.exponent() + edelta;
    let mb = b.mant();

    let d = (ea - eb) as u32;
    if d > 12 {
        // hypot(a, b) - a < b^2 / (2 * a) < a * 2^-25, which is less
        // than half an ULP of a.
        return a;
    }

    // hypot(a, b) = sqrt(Q) * 2^(ea - 23), with
    // Q = ma^2 + mb^2 * 2^(-2 * d)
    // 2^46 <= Q < 2^49
    //
    // p = 8 * Q, truncated to an integer, with the bits of
    // `mb^2 * 2^(-2 * d)` that are shifted out collapsed into a sticky bit
    // (they are non-zero iff 2 * tz(mb) + 3 < 2 * d). All the values that
    // `p` is compared with below are even integers, so the comparisons
    // give the same result as with the exact 8 * Q.
    let pa = (u64::from(ma) * u64::from(ma)) << 3;
    let pb = (u64::from(mb) * u64::from(mb)) << 3;
    let sticky = u64::from(2 * mb.trailing_zeros() + 3 < 2 * d);
    let p = pa + ((pb >> (2 * d)) | sticky);

    // m = Q * 2^-(16 + 2 * t), as 2.30 fixed point, so 1 <= m < 4
    let t = u32::from(p >= (1 << 51));
    let m = (p >> (19 + 2 * t)) as u32;

    // s ~= sqrt(Q) * 2^-t = sqrt(m) * 2^23, rounded to an integer
    // (2^23 <= s <= 2^24), with an error less than 1/8 before rounding
    let (_, s) = rsqrt_sqrt_32(m);
    let s = u64::from((s + (1 << 6)) >> 7);

    // `s` is the correctly rounded result or differs from it by one.
    // The correctly rounded result satisfies
    // (s - 1/2)^2 * 4^t <= Q <= (s + 1/2)^2 * 4^t (with ties to even),
    // that is, with k = 2 * t + 3,
    // -s * 2^k <= p - s^2 * 2^k - 2^(k - 2) <= s * 2^k.
    let k = 2 * t + 3;
    let rem = p.wrapping_sub((s * s) << k) as i64 - (1 << (k - 2));
    let lim = (s << k) as i64;
    let odd = (s & 1) as i64;
    let s = s + u64::from(rem + odd > lim) - u64::from(rem - odd < -lim);

    // hypot(a, b) = s * 2^(ea + t - 23), adding `s` (which has the
    // implicit bit) to the exponent minus one.
    let raw = (u64::from(f32::exp_to_raw_exp(ea + t as i16) - 1) << 23) + s;
    if raw >= u64::from(f32::INFINITY.to_bits()) {
        f32::INFINITY
    } else {
        f32::from_bits(raw as u32)
    }
}

/// Returns `sqrt(n)` rounded to the nearest integer, for `0 < n < 2^48`.
fn isqrt_round(n: u64) -> u32 {
    // q = n * 4^k, with 2^62 <= q < 2^64 (and k >= 8)
    let k = n.leading_zeros() / 2;
    let q = n << (2 * k);

    // s ~= sqrt(n) = sqrt(q) * 2^-k, rounded to an integer, with an error
    // less than 1/8 before rounding (since sqrt(n) < 2^24)
    let (_, s) = rsqrt_sqrt_32((q >> 32) as u32);
    let shift = k - 1;
    let s = (s + (1 << (shift - 1))) >> shift;

    // `s` is the correctly rounded result or differs from it by one.
    // The correctly rounded result satisfies (s - 1/2)^2 < n < (s + 1/2)^2,
    // which (since n and s are integers) is equivalent to -s < n - s^2 <= s.
    let rem = n as i64 - i64::from(s) * i64::from(s);
    s + u32::from(rem > i64::from(s)) - u32::from(rem <= -i64::from(s))
}
