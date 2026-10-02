use super::sqrt::sqrt_approx;
use crate::traits::Float as _;

impl crate::generic::Hypot for f64 {
    #[inline]
    fn hypot_finite(x: Self, y: Self) -> Self {
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
            return f64::from_bits(a);
        }

        if a < f64::MIN_POSITIVE.to_bits() {
            // Both are subnormal, so hypot(a, b) = sqrt(A^2 + B^2) * 2^-1074,
            // where A and B are the raw bits of `a` and `b`. The result is
            // less than 2^-1021, so its ULP is 2^-1074 and its raw bits are
            // sqrt(A^2 + B^2) rounded to an integer.
            return f64::from_bits(isqrt_round(
                u128::from(a) * u128::from(a) + u128::from(b) * u128::from(b),
            ));
        }

        // a = ma * 2^(ea - 52), with 2^52 <= ma < 2^53
        // b = mb * 2^(eb - 52), with 2^52 <= mb < 2^53
        let a = f64::from_bits(a);
        let ea = a.exponent();
        let ma = a.mant();
        let (b, edelta) = f64::from_bits(b).normalize_arg();
        let eb = b.exponent() + edelta;
        let mb = b.mant();

        let d = (ea - eb) as u32;
        if d > 27 {
            // hypot(a, b) - a < b^2 / (2 * a) < a * 2^-54, which is less
            // than half an ULP of a.
            return a;
        }

        // hypot(a, b) = sqrt(Q) * 2^(ea - 52), with
        // Q = ma^2 + mb^2 * 2^(-2 * d)
        // 2^104 <= Q < 2^107
        //
        // p = 8 * Q, truncated to an integer, with the bits of
        // `mb^2 * 2^(-2 * d)` that are shifted out collapsed into a sticky bit
        // (they are non-zero iff 2 * tz(mb) + 3 < 2 * d). All the values that
        // `p` is compared with below are even integers, so the comparisons
        // give the same result as with the exact 8 * Q.
        let pa = (u128::from(ma) * u128::from(ma)) << 3;
        let pb = (u128::from(mb) * u128::from(mb)) << 3;
        let sticky = u128::from(2 * mb.trailing_zeros() + 3 < 2 * d);
        let p = pa + ((pb >> (2 * d)) | sticky);

        // m = Q * 2^-(42 + 2 * t), as 2.62 fixed point, so 1 <= m < 4
        let t = u32::from(p >= (1 << 109));
        let m = (p >> (45 + 2 * t)) as u64;

        // s ~= sqrt(Q) * 2^-t = sqrt(m) * 2^52, rounded to an integer
        // (2^52 <= s <= 2^53), with an error less than 1/8 before rounding
        let s = (sqrt_approx(m) + (1 << 8)) >> 9;

        // `s` is the correctly rounded result or differs from it by one.
        // The correctly rounded result satisfies
        // (s - 1/2)^2 * 4^t <= Q <= (s + 1/2)^2 * 4^t (with ties to even),
        // that is, with k = 2 * t + 3,
        // -s * 2^k <= p - s^2 * 2^k - 2^(k - 2) <= s * 2^k.
        // `p - s^2 * 2^k` is small, so it can be calculated with wrapping
        // 64-bit arithmetic.
        let k = 2 * t + 3;
        let rem = (p as u64).wrapping_sub(s.wrapping_mul(s) << k) as i64 - (1 << (k - 2));
        let lim = (s << k) as i64;
        let odd = (s & 1) as i64;
        let s = s + u64::from(rem + odd > lim) - u64::from(rem - odd < -lim);

        // hypot(a, b) = s * 2^(ea + t - 52), adding `s` (which has the
        // implicit bit) to the exponent minus one.
        let raw = (u64::from(f64::exp_to_raw_exp(ea + t as i16) - 1) << 52) + s;
        if raw >= f64::INFINITY.to_bits() {
            f64::INFINITY
        } else {
            f64::from_bits(raw)
        }
    }
}

/// Returns `sqrt(n)` rounded to the nearest integer, for `0 < n < 2^106`.
fn isqrt_round(n: u128) -> u64 {
    // q = n * 4^k, with 2^126 <= q < 2^128 (and k >= 11)
    let k = n.leading_zeros() / 2;
    let q = n << (2 * k);

    // s ~= sqrt(n) = sqrt(q) * 2^-k, rounded to an integer, with an error
    // less than 1/8 before rounding (since sqrt(n) < 2^53)
    let shift = k - 2;
    let s = (sqrt_approx((q >> 64) as u64) + (1 << (shift - 1))) >> shift;

    // `s` is the correctly rounded result or differs from it by one.
    // The correctly rounded result satisfies (s - 1/2)^2 < n < (s + 1/2)^2,
    // which (since n and s are integers) is equivalent to -s < n - s^2 <= s.
    let rem = n.wrapping_sub(u128::from(s) * u128::from(s)) as i128;
    s + u64::from(rem > i128::from(s)) - u64::from(rem <= -i128::from(s))
}
