//! Trigonometric functions for `f32`.
//!
//! They use the same table and formulas as the `f64` functions (see
//! `crate::f64::trigonometric`), but the reduced argument is a single `f64`
//! and the evaluation uses plain `f64` arithmetic, which also allows a
//! cheaper tangent. The argument is reduced to `x = n * π/64 + b`, where `n`
//! is an integer (only `n mod 128` is used) and `|b| <= ~π/128`:
//!
//! * Radians: `|x| < 2^21` is reduced by `reduce_rad`, larger arguments use
//!   the `f64` reduction.
//! * Degrees and half revolutions: `|x| = n * 90/32 + t` and
//!   `|x| = (n + t)/64`, with `t` exact (`reduce_deg` and
//!   `reduce_half_revs`), followed by the conversion to radians. Huge
//!   arguments are first replaced with a smaller number with the same
//!   remainder modulo 360 (degrees) or with zero (half revolutions).
//!
//! With `a = n * π/64`, `s = sin(a)` and `c = cos(a)` (from the `f64` table):
//!
//! * `sin(a + b) = (s + c*b) + (s * (cos(b) - 1) + c * (sin(b) - b))`
//! * `cos(a + b) = (c - s*b) + (c * (cos(b) - 1) - s * (sin(b) - b))`
//! * `tan(a + b) = (s + c*tan(b)) / (c - s*tan(b))`
//!
//! All the functions reduce `|x|` (the odd ones flip the sign of the result
//! when `x` is negative), which gives exactly `+0`, `±1` or `±inf` (with the
//! right signs of zero) when `x` is a multiple of a right angle in degrees or
//! half revolutions.

use crate::f64::{FRAC_64_PI, reduce_rad_large_sum, round_u8, sin_cos_pi_64};
use crate::traits::Float;

// GENERATE: consts f64 FRAC_PI_180 PI
const FRAC_PI_180: f64 = f64::from_bits(0x3F91DF46A2529D39); // 1.7453292519943295e-2
const PI: f64 = f64::from_bits(0x400921FB54442D18); // 3.141592653589793e0

/// `π/64` (exact scaling of `PI`)
const FRAC_PI_64: f64 = PI * (1.0 / 64.0);

impl crate::generic::Trigonometric for f32 {
    #[inline]
    fn sin_finite(x: Self) -> Self {
        let (n, y) = reduce_rad(x);
        with_sign_of(sin_eval(n, y), x)
    }

    #[inline]
    fn cos_finite(x: Self) -> Self {
        let (n, y) = reduce_rad(x);
        cos_eval(n, y) as f32
    }

    #[inline]
    fn sin_cos_finite(x: Self) -> (Self, Self) {
        let (n, y) = reduce_rad(x);
        let (sin, cos) = sin_cos_eval(n, y);
        (with_sign_of(sin, x), cos as f32)
    }

    #[inline]
    fn tan_finite(x: Self) -> Self {
        let (n, y) = reduce_rad(x);
        with_sign_of(tan_eval(n, y), x)
    }

    #[inline]
    fn sind_finite(x: Self) -> Self {
        let (n, y) = reduce_deg(x);
        with_sign_of(sin_eval(n, y * FRAC_PI_180), x)
    }

    #[inline]
    fn cosd_finite(x: Self) -> Self {
        let (n, y) = reduce_deg(x);
        cos_eval(n, y * FRAC_PI_180) as f32
    }

    #[inline]
    fn sind_cosd_finite(x: Self) -> (Self, Self) {
        let (n, y) = reduce_deg(x);
        let (sin, cos) = sin_cos_eval(n, y * FRAC_PI_180);
        (with_sign_of(sin, x), cos as f32)
    }

    #[inline]
    fn tand_finite(x: Self) -> Self {
        let (n, y) = reduce_deg(x);
        with_sign_of(tan_eval(n, y * FRAC_PI_180), x)
    }

    #[inline]
    fn sinpi_finite(x: Self) -> Self {
        let (n, y) = reduce_half_revs(x);
        with_sign_of(sin_eval(n, y * FRAC_PI_64), x)
    }

    #[inline]
    fn cospi_finite(x: Self) -> Self {
        let (n, y) = reduce_half_revs(x);
        cos_eval(n, y * FRAC_PI_64) as f32
    }

    #[inline]
    fn sinpi_cospi_finite(x: Self) -> (Self, Self) {
        let (n, y) = reduce_half_revs(x);
        let (sin, cos) = sin_cos_eval(n, y * FRAC_PI_64);
        (with_sign_of(sin, x), cos as f32)
    }

    #[inline]
    fn tanpi_finite(x: Self) -> Self {
        let (n, y) = reduce_half_revs(x);
        with_sign_of(tan_eval(n, y * FRAC_PI_64), x)
    }
}

/// Rounds `v` to `f32` and flips its sign if `x` is negative.
#[inline]
fn with_sign_of(v: f64, x: f32) -> f32 {
    let v = v as f32;
    f32::from_bits(v.to_bits() ^ (x.to_bits() & <f32 as Float>::SIGN_MASK))
}

/// Reduces the angle argument `|x|` (in radians) to `(n, b)` such that:
/// * `|b| <= ~π/128`
/// * `|x| = 2*π*M + π/64*n + b`
/// * `M` is an integer
///
/// `n` is only meaningful modulo 128.
#[inline]
fn reduce_rad(x: f32) -> (u8, f64) {
    let xa = f64::from(x).abs();
    if x.exponent() < 21 {
        // GENERATE: split_const FRAC_PI_64_F PI -6 27 53
        const FRAC_PI_64_F0: f64 = f64::from_bits(0x3FA921FB54000000); // 4.9087385181337595e-2
        const FRAC_PI_64_F1: f64 = f64::from_bits(0x3DC10B4611A62633); // 3.100292436501689e-11

        // |x| < 2^21, so 0 <= n < 2^25.4 and, since `FRAC_PI_64_F0` has 27
        // significant bits, `nf * FRAC_PI_64_F0` is exact, and so is
        // `xa - nf * FRAC_PI_64_F0` (by Sterbenz lemma, or because `nf = 0`).
        // The error of the split (below `2^-88 * n`) and the rounding of
        // `nf * FRAC_PI_64_F1` (below `2^-87.9 * n`) add an absolute error
        // below `2^-61.5`, and the last subtraction a relative error below
        // `u`. When `n` is a nonzero multiple of 32, the distance from `x` to
        // `n * π/64` is at least 2^-24.8 when |x| >= 2^18, 2^-25.8 when
        // 2^15 <= |x| < 2^18 and 2^-27.8 when |x| < 2^15 (found by
        // exhaustive search), so the relative error of `b` is below 2^-38.8.
        let (nf, n) = round_u8(xa * FRAC_64_PI);
        let b = (xa - nf * FRAC_PI_64_F0) - nf * FRAC_PI_64_F1;
        (n, b)
    } else {
        reduce_rad_large_sum(xa)
    }
}

/// Reduces the angle argument `|x|` (in degrees) to `(n, t)` such that:
/// * `|t| <= ~45/32`
/// * `|x| = 360*M + 90/32*n + t`
/// * `M` is an integer
/// * `t` is exact
///
/// `n` is only meaningful modulo 128.
#[inline]
fn reduce_deg(x: f32) -> (u8, f64) {
    // `STEP` has 6 significant bits.
    const STEP: f64 = 90.0 / 32.0;

    let y = if x.exponent() < 26 {
        f64::from(x).abs()
    } else {
        deg_mod_360(x)
    };

    // 0 <= y < 2^33, so 0 <= n < 2^32 and `nf * STEP` is exact.
    // `y - nf * STEP` is exact by Sterbenz lemma when n >= 2 or y >= 2.
    // Otherwise, both operands are multiples of `ulp(y)` and the result is
    // less than 2 in magnitude (or `nf = 0`).
    let (nf, n) = round_u8(y * (1.0 / STEP));
    (n, y - nf * STEP)
}

/// Returns `y` such that `0 <= y < 2^33` and `y = |x| (mod 360)`, for
/// `|x| >= 2^26`.
///
/// Kept out of line: when inlined, the `f32` to `f64` conversion in the fast
/// path of `reduce_deg` tends to be allocated to a register other than its
/// source, and SSE `cvtss2sd` then depends on the stale value of that
/// register (e.g., the result of a previous call), serializing calls.
#[cold]
fn deg_mod_360(x: f32) -> f64 {
    // |x| = m * 2^s, with `m` a 24-bit integer and s >= 3, so `|x| mod 360`
    // is the remainder of `m * (2^s mod 360)`, with
    // 2^s mod 360 = 8 * (2^(s - 3) mod 45) and 2^12 = 1 (mod 45).
    // m * (2^s mod 360) < 2^24 * 360 < 2^33 is exact.
    let m = x.mant();
    let s = (x.exponent() - 23) as u32;
    let pow2 = 8 * ((1u32 << ((s - 3) % 12)) % 45);
    f64::from(m) * f64::from(pow2)
}

/// Reduces the angle argument `|x|` (in half revolutions) to `(n, t)` such
/// that:
/// * `|t| <= 1/2`
/// * `|x| = 2*M + (n + t)/64`
/// * `M` is an integer
/// * `t` is exact
///
/// `n` is only meaningful modulo 128.
#[inline]
fn reduce_half_revs(x: f32) -> (u8, f64) {
    // |x| >= 2^24 is an even integer, so it can be replaced with zero.
    let y = if x.exponent() < 24 {
        f64::from(x).abs()
    } else {
        0.0
    };

    // 0 <= y < 2^24, so `u = y * 64` is exact and less than 2^30, and so is
    // `u - nf`.
    let u = y * 64.0;
    let (nf, n) = round_u8(u);
    (n, u - nf)
}

/// Returns `sin(π * x)` without rounding it to `f32` (used by the reflection
/// formula of `ln_gamma`).
#[inline]
pub(crate) fn sinpi_f64(x: f32) -> f64 {
    let (n, t) = reduce_half_revs(x);
    let s = sin_eval(n, t * FRAC_PI_64);
    if x.is_sign_negative() { -s } else { s }
}

/// Returns `(cos(b) - 1, sin(b) - b)`, `|b| <= 0.0247`.
#[inline]
fn sin_cos_poly(b: f64) -> (f64, f64) {
    // GENERATE: sin_poly f64 2 0.0247
    const K3: f64 = f64::from_bits(0xBFC555555543D805); // -1.6666666663485344e-1
    const K5: f64 = f64::from_bits(0x3F8110FA72342413); // 8.333164797038639e-3

    // GENERATE: cos_poly f64 1 0.0247
    const K4: f64 = f64::from_bits(0x3FA5553BE97F5E6D); // 4.166590905486776e-2

    let b2 = b * b;
    let cm1 = b2 * (-0.5 + b2 * K4);
    let smb = (b * b2) * (K3 + b2 * K5);
    (cm1, smb)
}

/// Returns `sin(n * π/64 + b)`.
#[inline]
fn sin_eval(n: u8, b: f64) -> f64 {
    let (s, c) = sin_cos_pi_64(n);
    let (cm1, smb) = sin_cos_poly(b);
    (s + c * b) + (s * cm1 + c * smb)
}

/// Returns `cos(n * π/64 + b)`.
#[inline]
fn cos_eval(n: u8, b: f64) -> f64 {
    let (s, c) = sin_cos_pi_64(n);
    let (cm1, smb) = sin_cos_poly(b);
    (c - s * b) + (c * cm1 - s * smb)
}

/// Returns `(sin(n * π/64 + b), cos(n * π/64 + b))`.
///
/// Evaluated as `(s * cos(b) + c * sin(b), c * cos(b) - s * sin(b))`, sharing
/// `cos(b)` and `sin(b)`.
#[inline]
fn sin_cos_eval(n: u8, b: f64) -> (f64, f64) {
    let (s, c) = sin_cos_pi_64(n);
    let (cm1, smb) = sin_cos_poly(b);
    let cb = 1.0 + cm1;
    let sb = b + smb;
    (s * cb + c * sb, c * cb - s * sb)
}

/// Returns `tan(n * π/64 + b)`.
#[inline]
fn tan_eval(n: u8, b: f64) -> f64 {
    // GENERATE: tan_poly f64 2 0.0247
    const K3: f64 = f64::from_bits(0x3FD555554C09CA81); // 3.3333332467660887e-1
    const K5: f64 = f64::from_bits(0x3FC11291BAE7B019); // 1.333791887876721e-1

    let (s, c) = sin_cos_pi_64(n);
    // t ~= tan(b)
    let b2 = b * b;
    let t = b + (b * b2) * (K3 + b2 * K5);
    // tan(a + b) = (sin(a) + cos(a) * tan(b)) / (cos(a) - sin(a) * tan(b))
    (s + c * t) / (c - s * t)
}
