//! Trigonometric functions for `f64`, and the parts that `f32` shares.
//!
//! The argument is reduced to `x = n * π/64 + b` (or the equivalent in
//! degrees or half revolutions), where `n` is an integer (only `n mod 128` is
//! used) and `|b| <= ~π/128`. With `a = n * π/64`, `sin(a)` and
//! `cos(a) = sin(a + π/2)` are taken from `SIN_PI_TBL`, and
//!
//! `sin(a + b) = sin(a) + cos(a) * b + sin(a) * (cos(b) - 1) + cos(a) * (sin(b) - b)`
//!
//! where `cos(b) - 1` and `sin(b) - b` are approximated with short polynomials.
//! The tangent is the quotient of the sine and the cosine.
//!
//! Reductions:
//! * Radians: Cody-Waite with splits of π/64 (`reduce_rad`), or Payne-Hanek
//!   for huge arguments (`reduce_rad_huge`). `b` is a double-word.
//! * Degrees and half revolutions: `x = n * 90/32 + t` and `x = (n + t)/64`,
//!   with `t` exact (`reduce_deg` and `reduce_half_revs`), then `t` is
//!   converted to radians. Huge arguments are first reduced modulo 360 and 2
//!   with integer arithmetic.
//!
//! The odd functions in degrees and half revolutions reduce `|x|` and flip the
//! sign of the result when `x` is negative. When `x` is a multiple of a right
//! angle (`t = 0` and `n` multiple of 32), the evaluation yields exactly `+0`,
//! `±1` or `±inf` with the signs of zero required by `sin(|x|)`, `cos(|x|)`
//! and `tan(|x|)`.

use super::f64x2::F64x2;
use crate::traits::Float as _;

// GENERATE: consts F64x2 FRAC_PI_180 PI
const FRAC_PI_180: F64x2 = F64x2::from_bits(0x3F91DF46A2529D39, 0x3C15C1D8BECDD291); // 1.745329251994329576923690768489e-2
const PI: F64x2 = F64x2::from_bits(0x400921FB54442D18, 0x3CA1A62633145C07); // 3.141592653589793238462643383280e0

// GENERATE: consts f64 FRAC_2_PI
const FRAC_2_PI: f64 = f64::from_bits(0x3FE45F306DC9C883); // 6.366197723675814e-1

/// `64 / π`
pub(crate) const FRAC_64_PI: f64 = FRAC_2_PI * 32.0;

impl crate::generic::Trigonometric for f64 {
    #[inline]
    fn sin_finite(x: Self) -> Self {
        let (n, y) = reduce_rad(x);
        y.sin_k(n).to_f64()
    }

    #[inline]
    fn cos_finite(x: Self) -> Self {
        let (n, y) = reduce_rad(x);
        y.cos_k(n).to_f64()
    }

    #[inline]
    fn sin_cos_finite(x: Self) -> (Self, Self) {
        let (n, y) = reduce_rad(x);
        (y.sin_k(n).to_f64(), y.cos_k(n).to_f64())
    }

    #[inline]
    fn tan_finite(x: Self) -> Self {
        if x.exponent() <= -40 {
            // Tiny, tan(x) ~= x
            x
        } else {
            let (n, y) = reduce_rad(x);
            y.tan_k(n)
        }
    }

    #[inline]
    fn sind_finite(x: Self) -> Self {
        if x.exponent() < -960 {
            // Avoid subnormals in `F64x2`, which break accuracy.
            (F64x2::new1(x * f64::exp2i_fast(106)) * FRAC_PI_180).scalbn_to_f64(-106)
        } else {
            let (n, y) = reduce_deg(x.abs());
            with_sign_of(deg_to_rad(y).sin_k(n).to_f64(), x)
        }
    }

    #[inline]
    fn cosd_finite(x: Self) -> Self {
        let (n, y) = reduce_deg(x.abs());
        deg_to_rad(y).cos_k(n).to_f64()
    }

    #[inline]
    fn sind_cosd_finite(x: Self) -> (Self, Self) {
        if x.exponent() < -960 {
            // Avoid subnormals in `F64x2`, which break accuracy.
            let sin = (F64x2::new1(x * f64::exp2i_fast(106)) * FRAC_PI_180).scalbn_to_f64(-106);
            (sin, 1.0)
        } else {
            let (n, y) = reduce_deg(x.abs());
            let z = deg_to_rad(y);
            (with_sign_of(z.sin_k(n).to_f64(), x), z.cos_k(n).to_f64())
        }
    }

    #[inline]
    fn tand_finite(x: Self) -> Self {
        if x.exponent() <= -40 {
            // Avoid subnormals in `F64x2`, which break accuracy.
            (F64x2::new1(x * f64::exp2i_fast(106)) * FRAC_PI_180).scalbn_to_f64(-106)
        } else {
            let (n, y) = reduce_deg(x.abs());
            with_sign_of(deg_to_rad(y).tan_k(n), x)
        }
    }

    #[inline]
    fn sinpi_finite(x: Self) -> Self {
        if x.exponent() < -960 {
            // Avoid subnormals in `F64x2`, which break accuracy.
            (F64x2::new1(x * f64::exp2i_fast(106)) * PI).scalbn_to_f64(-106)
        } else {
            let (n, y) = reduce_half_revs(x.abs());
            with_sign_of(half_revs_to_rad(y).sin_k(n).to_f64(), x)
        }
    }

    #[inline]
    fn cospi_finite(x: Self) -> Self {
        let (n, y) = reduce_half_revs(x.abs());
        half_revs_to_rad(y).cos_k(n).to_f64()
    }

    #[inline]
    fn sinpi_cospi_finite(x: Self) -> (Self, Self) {
        if x.exponent() < -960 {
            // Avoid subnormals in `F64x2`, which break accuracy.
            let sin = (F64x2::new1(x * f64::exp2i_fast(106)) * PI).scalbn_to_f64(-106);
            (sin, 1.0)
        } else {
            let (n, y) = reduce_half_revs(x.abs());
            let z = half_revs_to_rad(y);
            (with_sign_of(z.sin_k(n).to_f64(), x), z.cos_k(n).to_f64())
        }
    }

    #[inline]
    fn tanpi_finite(x: Self) -> Self {
        if x.exponent() <= -40 {
            // Avoid subnormals in `F64x2`, which break accuracy.
            (F64x2::new1(x * f64::exp2i_fast(106)) * PI).scalbn_to_f64(-106)
        } else {
            let (n, y) = reduce_half_revs(x.abs());
            with_sign_of(half_revs_to_rad(y).tan_k(n), x)
        }
    }
}

/// Flips the sign of `v` if `x` is negative.
#[inline]
fn with_sign_of(v: f64, x: f64) -> f64 {
    f64::from_bits(v.to_bits() ^ (x.to_bits() & (1 << 63)))
}

/// Rounds `x` to the nearest integer (ties to even), returning it as
/// `f64` and its lowest 8 bits (in two's complement).
///
/// `|x|` must be less than 2^51.
#[inline]
pub(crate) fn round_i(x: f64) -> (f64, u8) {
    // 1.5 * 2^52, adding it rounds to an integer, keeping the integer in the
    // lowest bits of the mantissa.
    const MAGIC: f64 = 6755399441055744.0;

    let t = (x + MAGIC).purify();
    let r = (t - MAGIC).purify();
    (r, t.to_bits() as u8)
}

/// Returns `(sin(n * π/64), cos(n * π/64))`, rounded to `f64`.
#[inline]
pub(crate) fn sin_cos_pi_64(n: u8) -> (f64, f64) {
    (
        SIN_PI_TBL[usize::from(n & 127)].hi(),
        SIN_PI_TBL[usize::from(n.wrapping_add(32) & 127)].hi(),
    )
}

/// Reduces the angle argument `x` (in radians) to `(n, b)` such that:
/// * `|b| <= ~π/128`
/// * `x = 2*π*M + π/64*n + b`
/// * `M` is an integer
///
/// `n` is only meaningful modulo 128. The reduction of `-x` is `(-n, -b)`.
#[inline]
fn reduce_rad(x: f64) -> (u8, Reduced) {
    let x = x.purify();
    let absx = x.abs();
    if absx < 1.5 {
        // GENERATE: reduce_rad::split_pi FRAC_PI_64_S -6 48 53
        const FRAC_PI_64_S0: f64 = f64::from_bits(0x3FA921FB54442D00); // 4.908738521234035e-2
        const FRAC_PI_64_S1: f64 = f64::from_bits(0x3CA8469898CC5170); // 1.6844696431744122e-16

        // |n| <= 31 and `FRAC_PI_64_S0` has 48 significant bits, so
        // `nf * FRAC_PI_64_S0` is exact. `x - nf * FRAC_PI_64_S0` is exact:
        // * For |n| >= 2 or |x| >= 2^-5, by Sterbenz lemma.
        // * For |n| = 1 and |x| < 2^-5, both operands are multiples of
        //   `ulp(x)` and the result is less than 2^-5 <= 2^53 * ulp(x).
        // * For n = 0, trivially.
        //
        // `bl` is not small compared to `bh` (up to ~2^-47.4), but its
        // first and second order terms are taken into account by `sin_a`.
        let (nf, n) = round_i(x * FRAC_64_PI);

        let bh = (x - nf * FRAC_PI_64_S0).purify();
        let bl = -(nf * FRAC_PI_64_S1);

        (n, Reduced { bh, bl })
    } else if absx < 256.0 {
        // GENERATE: reduce_rad::split_pi FRAC_PI_64_M -6 40 40 53
        const FRAC_PI_64_M0: f64 = f64::from_bits(0x3FA921FB54442000); // 4.908738521231726e-2
        const FRAC_PI_64_M1: f64 = f64::from_bits(0x3D1A308D31318000); // 2.3261085876500796e-14
        const FRAC_PI_64_M2: f64 = f64::from_bits(0x3A98A2E03707344A); // 1.9900991135971718e-26

        // 31 <= |n| < 2^13, so `nf * FRAC_PI_64_M0` and `nf * FRAC_PI_64_M1`
        // are exact, and `x - nf * FRAC_PI_64_M0` is exact by Sterbenz lemma.
        let (nf, n) = round_i(x * FRAC_64_PI);

        let r0 = (x - nf * FRAC_PI_64_M0).purify();
        let r1 = F64x2::sub11(r0, (nf * FRAC_PI_64_M1).purify());
        let bl = r1.lo() - nf * FRAC_PI_64_M2;

        (n, Reduced { bh: r1.hi(), bl })
    } else {
        reduce_rad_large(x)
    }
}

/// Like `reduce_rad`, for `|x| >= 256`, returning `(n, b)`.
fn reduce_rad_large(x: f64) -> (u8, Reduced) {
    if x.exponent() < 39 {
        // GENERATE: reduce_rad::split_pi FRAC_PI_64_P -6 53 53 53
        const FRAC_PI_64_P0: f64 = f64::from_bits(0x3FA921FB54442D18); // 4.908738521234052e-2
        const FRAC_PI_64_P1: f64 = f64::from_bits(0x3C41A62633145C06); // 1.913510623667739e-18
        const FRAC_PI_64_P2: f64 = f64::from_bits(0x38FC1CD129024E09); // 3.3839271060059813e-34

        // 2^12 < |n| < 2^44
        let (nf, n) = round_i(x * FRAC_64_PI);

        // n * π/64 ~= (ph + pl) + (qh + ql) + nf * FRAC_PI_64_P2
        let p = F64x2::mul11(nf, FRAC_PI_64_P0);
        let q = F64x2::mul11(nf, FRAC_PI_64_P1);

        // Exact by Sterbenz lemma
        let rh = (x - p.hi()).purify();
        // `rh` is zero or a multiple of `ulp(p.hi())/2`, so `|rh| >= |p.lo()|`
        // and `fast_add11` is exact.
        let r = F64x2::fast_add11(rh, -p.lo());

        let y = F64x2::sub11(r.hi(), q.hi());
        let bl = y.lo() + (r.lo() - (q.lo() + nf * FRAC_PI_64_P2));

        (n, Reduced { bh: y.hi(), bl })
    } else {
        reduce_rad_huge(x)
    }
}

/// Like `reduce_rad_large`, with `b` rounded to a single `f64` (used by
/// `f32`).
#[inline]
pub(crate) fn reduce_rad_large_sum(x: f64) -> (u8, f64) {
    let (n, b) = reduce_rad_large(x);
    (n, b.bh + b.bl)
}

/// Like `reduce_rad_large`, for `|x| >= 2^39`.
///
/// With `|x| = m * 2^e` (`m` a 53-bit integer), `x * 64/π` is computed modulo
/// 128 as `m * G * 2^7`, where `G = frac(2^(e-1) / π)` is taken from the bits of
/// `1/π` starting at position `e` (the discarded leading bits contribute
/// multiples of 128). 192 bits of `G` are used, so the truncation error of
/// `x * 64/π` is less than `2^53 * 2^-185 = 2^-132`, and the fraction is
/// kept with 128 bits (error less than 2^-128). Since the reduced argument is at
/// least ~2^-61 when it is close to a multiple of π/2, this is accurate enough.
#[cold]
fn reduce_rad_huge(x: f64) -> (u8, Reduced) {
    // GENERATE: reduce_rad::frac_1_pi_bits 63 20
    // 2^-63 / π ~= sum(FRAC_1_PI_BITS[i] * 2^(-64 * (i + 1)))
    static FRAC_1_PI_BITS: [u64; 20] = [
        0x0000000000000000,
        0xA2F9836E4E441529,
        0xFC2757D1F534DDC0,
        0xDB6295993C439041,
        0xFE5163ABDEBBC561,
        0xB7246E3A424DD2E0,
        0x06492EEA09D1921C,
        0xFE1DEB1CB129A73E,
        0xE88235F52EBB4484,
        0xE99C7026B45F7E41,
        0x3991D639835339F4,
        0x9C845F8BBDF9283B,
        0x1FF897FFDE05980F,
        0xEF2F118B5A0A6D1F,
        0x6D367ECF27CB09B7,
        0x4F463F669E5FEA2D,
        0x7527BAC7EBE5F17B,
        0x3D0739F78A5292EA,
        0x6BFB5FB11F8D5D08,
        0x56033046FC7B6BAB,
    ];

    let m = u128::from(x.mant());
    // |x| = m * 2^e. With `1/π = sum(b[i] * 2^-i)`, the first bit of `G` is
    // `b[e]`, which is bit `e + 62` (counting from 0, starting at the most
    // significant bit of the first word) of `FRAC_1_PI_BITS`.
    let pos = (i32::from(x.exponent()) - 52 + 62) as usize;
    let (w, s) = (pos / 64, pos % 64);
    let window = |i: usize| -> u64 {
        let pair = (u128::from(FRAC_1_PI_BITS[i]) << 64) | u128::from(FRAC_1_PI_BITS[i + 1]);
        ((pair << s) >> 64) as u64
    };

    // q = m * G (192 bits of G), only the lowest 192 bits of the result are
    // needed: bits 185..192 are the integer part (modulo 128) of x * 64/π,
    // and bits 0..185 the fraction.
    let t = m * u128::from(window(w + 2));
    let q0 = t as u64;
    let t = m * u128::from(window(w + 1)) + (t >> 64);
    let q1 = t as u64;
    let t = m * u128::from(window(w)) + (t >> 64);
    let q2 = t as u64;

    let mut n = (q2 >> 57) as u8;
    // Highest 128 bits of the fraction
    let mut frac =
        (u128::from(q2 & ((1 << 57) - 1)) << 71) | (u128::from(q1) << 7) | u128::from(q0 >> 57);
    let neg = (frac >> 127) != 0;
    if neg {
        // Round to nearest, the fraction becomes negative
        n = n.wrapping_add(1);
        frac = frac.wrapping_neg();
    }

    // |fraction| = frac * 2^-128. Convert it to a non-normalized
    // double-double (fh, fl), dropping at most 22 bits (relative error
    // < 2^-105). When `n` is a multiple of 32, |x - n * π/64| >= 2^-60.89
    // (the minimum for all finite `f64`), so `frac >= 2^71`. Otherwise the
    // result is not small, so a small absolute error is enough (the masking
    // of `lz` avoids an overflowing shift in the impossible case
    // `frac == 0`).
    let lz = frac.leading_zeros() & 127;
    let frac = frac << lz;
    let fh = ((frac >> 75) as u64) as f64 * f64::exp2i_fast(-53 - lz as i16);
    let fl = (((frac >> 22) as u64) & ((1 << 53) - 1)) as f64 * f64::exp2i_fast(-106 - lz as i16);

    // b = fraction * π/64
    let p = F64x2::mul11(fh, PI.hi() * (1.0 / 64.0));
    let bl = p.lo() + (fh * (PI.lo() * (1.0 / 64.0)) + fl * (PI.hi() * (1.0 / 64.0)));
    let (bh, bl) = if neg { (-p.hi(), -bl) } else { (p.hi(), bl) };

    if x.is_sign_negative() {
        (n.wrapping_neg(), Reduced { bh: -bh, bl: -bl })
    } else {
        (n, Reduced { bh, bl })
    }
}

/// Reduces the angle argument `x` (in degrees) to `(n, t)` such that:
/// * `|t| <= ~45/32`
/// * `x = 360*M + 90/32*n + t` (`|x|` instead of `x` when `|x| >= 2^45`)
/// * `M` is an integer
/// * `t` is exact
///
/// `n` is only meaningful modulo 128. When `|x| < 2^45`, the reduction of
/// `-x` is `(-n, -t)`.
#[inline]
fn reduce_deg(x: f64) -> (u8, f64) {
    let x = x.purify();
    let y = if x.exponent() < 45 { x } else { deg_mod_360(x) };
    reduce_deg_small(y)
}

/// `reduce_deg` for `|y| < 2^45`.
#[inline]
fn reduce_deg_small(y: f64) -> (u8, f64) {
    // `STEP` has 6 significant bits and |n| < 2^44, so `nf * STEP` is exact.
    // `y - nf * STEP` is exact by Sterbenz lemma when |n| >= 2 or |y| >= 2.
    // Otherwise, both operands are multiples of `ulp(y)` and the result is
    // less than 2 in magnitude (or `nf = 0`).
    const STEP: f64 = 90.0 / 32.0;

    let (nf, n) = round_i(y * (1.0 / STEP));
    (n, y - nf * STEP)
}

/// Returns `|x| mod 360` (exactly), for `|x| >= 2^45`.
#[cold]
fn deg_mod_360(x: f64) -> f64 {
    // |x| = m * 2^s, with `m` a 53-bit integer and s >= -7, so
    // `X = |x| * 2^7 = m * 2^p` (p = s + 7) is an integer and
    // |x| mod 360 = (X mod 46080) / 2^7, with 46080 = 360 * 2^7 = 2^10 * 45.
    // 2^p mod 46080 is 2^p when p < 10, and 2^10 * (2^(p - 10) mod 45)
    // otherwise, where 2^12 = 1 (mod 45).
    let m = x.mant();
    let p = (x.exponent() - 45) as u32;
    let pow2 = if p < 10 {
        1 << p
    } else {
        (1 << 10) * ((1u64 << ((p - 10) % 12)) % 45)
    };
    let r = (m % 46080) * pow2 % 46080;
    r as f64 * (1.0 / 128.0)
}

/// Converts `t` from degrees to radians.
#[inline]
fn deg_to_rad(t: f64) -> Reduced {
    let p = F64x2::mul11(t, FRAC_PI_180.hi());
    let bl = p.lo() + t * FRAC_PI_180.lo();
    Reduced { bh: p.hi(), bl }
}

/// Reduces the angle argument `x` (in half revolutions) to `(n, t)` such
/// that:
/// * `|t| <= 1/2`
/// * `x = 2*M + (n + t)/64` (`|x|` instead of `x` when `|x| >= 2^45`)
/// * `M` is an integer
/// * `t` is exact
///
/// `n` is only meaningful modulo 128. When `|x| < 2^45`, the reduction of
/// `-x` is `(-n, -t)`.
#[inline]
fn reduce_half_revs(x: f64) -> (u8, f64) {
    let x = x.purify();
    let y = if x.exponent() < 45 {
        x
    } else {
        half_revs_mod_2(x)
    };
    reduce_half_revs_small(y)
}

/// `reduce_half_revs` for `|y| < 2^45`.
#[inline]
fn reduce_half_revs_small(y: f64) -> (u8, f64) {
    // `u = y * 64` is exact, and so is `u - nf`.
    let u = y * 64.0;
    let (nf, n) = round_i(u);
    (n, u - nf)
}

/// Returns `|x| mod 2` (exactly), for `|x| >= 2^45`.
#[cold]
fn half_revs_mod_2(x: f64) -> f64 {
    // |x| = m * 2^s, with `m` a 53-bit integer and s >= -7, so
    // `X = |x| * 2^7 = m * 2^p` (p = s + 7) is an integer and
    // |x| mod 2 = (X mod 256) / 2^7, which is zero when p >= 8.
    let p = (x.exponent() - 45) as u32;
    let r = if p < 8 { (x.mant() << p) & 255 } else { 0 };
    r as f64 * (1.0 / 128.0)
}

/// Converts `t` from units of 1/64 half revolution to radians.
#[inline]
fn half_revs_to_rad(t: f64) -> Reduced {
    let p = F64x2::mul11(t, PI.hi() * (1.0 / 64.0));
    let bl = p.lo() + t * (PI.lo() * (1.0 / 64.0));
    Reduced { bh: p.hi(), bl }
}

/// Returns `sin(π * x)` as a double-word, for finite `x` with
/// `|x| >= 2^-960` (used by the reflection formula of `ln_gamma`).
#[inline]
pub(crate) fn sinpi_f64x2(x: f64) -> F64x2 {
    let (n, t) = reduce_half_revs(x.abs());
    let s = half_revs_to_rad(t).sin_k(n).normalize();
    if x.is_sign_negative() { -s } else { s }
}

/// A reduced argument `b = bh + bl`, with `|b| <= ~π/128`.
///
/// `bh + bl` is not necessarily normalized, but `|bl| <= ~2^-47`.
struct Reduced {
    bh: f64,
    bl: f64,
}

// GENERATE: sin_pi_table F64x2 64
// SIN_PI_TBL[i] = sin(i * π / 64)
static SIN_PI_TBL: [F64x2; 128] = [
    F64x2::from_bits(0x0000000000000000, 0x0000000000000000), // 0
    F64x2::from_bits(0x3FA91F65F10DD814, 0xBC2912BD0D569A90), // 4.906767432741801425495497694268e-2
    F64x2::from_bits(0x3FB917A6BC29B42C, 0xBC3E2718D26ED688), // 9.801714032956060199419556388864e-2
    F64x2::from_bits(0x3FC2C8106E8E613A, 0x3C513000A89A11E0), // 1.467304744553617516588501296467e-1
    F64x2::from_bits(0x3FC8F8B83C69A60B, 0xBC626D19B9FF8D82), // 1.950903220161282678482848684770e-1
    F64x2::from_bits(0x3FCF19F97B215F1B, 0xBC642DEEF11DA2C4), // 2.429801799032638899482741620775e-1
    F64x2::from_bits(0x3FD294062ED59F06, 0xBC75D28DA2C4612D), // 2.902846772544623676361923758174e-1
    F64x2::from_bits(0x3FD58F9A75AB1FDD, 0xBC1EFDC0D58CF620), // 3.368898533922200506892532126191e-1
    F64x2::from_bits(0x3FD87DE2A6AEA963, 0xBC672CEDD3D5A610), // 3.826834323650897717284599840304e-1
    F64x2::from_bits(0x3FDB5D1009E15CC0, 0x3C65B362CB974183), // 4.275550934302820943209668568888e-1
    F64x2::from_bits(0x3FDE2B5D3806F63B, 0x3C5E0D891D3C6841), // 4.713967368259976485563876259053e-1
    F64x2::from_bits(0x3FE073879922FFEE, 0xBC8A5A014347406C), // 5.141027441932217265936938389688e-1
    F64x2::from_bits(0x3FE1C73B39AE68C8, 0x3C8B25DD267F6600), // 5.555702330196022247428308139485e-1
    F64x2::from_bits(0x3FE30FF7FCE17035, 0xBC6EFCC626F74A6F), // 5.956993044924333434670365288300e-1
    F64x2::from_bits(0x3FE44CF325091DD6, 0x3C68076A2CFDC6B3), // 6.343932841636454982151716132255e-1
    F64x2::from_bits(0x3FE57D69348CECA0, 0xBC875720992BFBB2), // 6.715589548470184006253768504274e-1
    F64x2::from_bits(0x3FE6A09E667F3BCD, 0xBC8BDD3413B26456), // 7.071067811865475244008443621048e-1
    F64x2::from_bits(0x3FE7B5DF226AAFAF, 0xBC70F537ACDF0AD7), // 7.409511253549590911756168974952e-1
    F64x2::from_bits(0x3FE8BC806B151741, 0xBC82C5E12ED1336D), // 7.730104533627369608109066097585e-1
    F64x2::from_bits(0x3FE9B3E047F38741, 0xBC830EE286712474), // 8.032075314806449098066765129631e-1
    F64x2::from_bits(0x3FEA9B66290EA1A3, 0x3C39F630E8B6DAC8), // 8.314696123025452370787883776179e-1
    F64x2::from_bits(0x3FEB728345196E3E, 0xBC8BC69F324E6D61), // 8.577286100002720699022699842848e-1
    F64x2::from_bits(0x3FEC38B2F180BDB1, 0xBC76E0B1757C8D07), // 8.819212643483550297127568636604e-1
    F64x2::from_bits(0x3FECED7AF43CC773, 0xBC5E7B6BB5AB58AE), // 9.039892931234433315862002972305e-1
    F64x2::from_bits(0x3FED906BCF328D46, 0x3C7457E610231AC2), // 9.238795325112867561281831893968e-1
    F64x2::from_bits(0x3FEE212104F686E5, 0xBC8014C76C126527), // 9.415440651830207784125094025995e-1
    F64x2::from_bits(0x3FEE9F4156C62DDA, 0x3C8760B1E2E3F81E), // 9.569403357322088649357978869803e-1
    F64x2::from_bits(0x3FEF0A7EFB9230D7, 0x3C752C7ADC6B4989), // 9.700312531945439926039842072861e-1
    F64x2::from_bits(0x3FEF6297CFF75CB0, 0x3C7562172A361FD3), // 9.807852804032304491261822361342e-1
    F64x2::from_bits(0x3FEFA7557F08A517, 0xBC87A0A8CA13571F), // 9.891765099647809734516737380162e-1
    F64x2::from_bits(0x3FEFD88DA3D12526, 0xBC887DF6378811C7), // 9.951847266721968862448369531095e-1
    F64x2::from_bits(0x3FEFF621E3796D7E, 0xBC6C57BC2E24AA15), // 9.987954562051723927147716047591e-1
    F64x2::from_bits(0x3FF0000000000000, 0x0000000000000000), // 1.000000000000000000000000000000e0
    F64x2::from_bits(0x3FEFF621E3796D7E, 0xBC6C57BC2E24AA15), // 9.987954562051723927147716047591e-1
    F64x2::from_bits(0x3FEFD88DA3D12526, 0xBC887DF6378811C7), // 9.951847266721968862448369531095e-1
    F64x2::from_bits(0x3FEFA7557F08A517, 0xBC87A0A8CA13571F), // 9.891765099647809734516737380162e-1
    F64x2::from_bits(0x3FEF6297CFF75CB0, 0x3C7562172A361FD3), // 9.807852804032304491261822361342e-1
    F64x2::from_bits(0x3FEF0A7EFB9230D7, 0x3C752C7ADC6B4989), // 9.700312531945439926039842072861e-1
    F64x2::from_bits(0x3FEE9F4156C62DDA, 0x3C8760B1E2E3F81E), // 9.569403357322088649357978869803e-1
    F64x2::from_bits(0x3FEE212104F686E5, 0xBC8014C76C126527), // 9.415440651830207784125094025995e-1
    F64x2::from_bits(0x3FED906BCF328D46, 0x3C7457E610231AC2), // 9.238795325112867561281831893968e-1
    F64x2::from_bits(0x3FECED7AF43CC773, 0xBC5E7B6BB5AB58AE), // 9.039892931234433315862002972305e-1
    F64x2::from_bits(0x3FEC38B2F180BDB1, 0xBC76E0B1757C8D07), // 8.819212643483550297127568636604e-1
    F64x2::from_bits(0x3FEB728345196E3E, 0xBC8BC69F324E6D61), // 8.577286100002720699022699842848e-1
    F64x2::from_bits(0x3FEA9B66290EA1A3, 0x3C39F630E8B6DAC8), // 8.314696123025452370787883776179e-1
    F64x2::from_bits(0x3FE9B3E047F38741, 0xBC830EE286712474), // 8.032075314806449098066765129631e-1
    F64x2::from_bits(0x3FE8BC806B151741, 0xBC82C5E12ED1336D), // 7.730104533627369608109066097585e-1
    F64x2::from_bits(0x3FE7B5DF226AAFAF, 0xBC70F537ACDF0AD7), // 7.409511253549590911756168974952e-1
    F64x2::from_bits(0x3FE6A09E667F3BCD, 0xBC8BDD3413B26456), // 7.071067811865475244008443621048e-1
    F64x2::from_bits(0x3FE57D69348CECA0, 0xBC875720992BFBB2), // 6.715589548470184006253768504274e-1
    F64x2::from_bits(0x3FE44CF325091DD6, 0x3C68076A2CFDC6B3), // 6.343932841636454982151716132255e-1
    F64x2::from_bits(0x3FE30FF7FCE17035, 0xBC6EFCC626F74A6F), // 5.956993044924333434670365288300e-1
    F64x2::from_bits(0x3FE1C73B39AE68C8, 0x3C8B25DD267F6600), // 5.555702330196022247428308139485e-1
    F64x2::from_bits(0x3FE073879922FFEE, 0xBC8A5A014347406C), // 5.141027441932217265936938389688e-1
    F64x2::from_bits(0x3FDE2B5D3806F63B, 0x3C5E0D891D3C6841), // 4.713967368259976485563876259053e-1
    F64x2::from_bits(0x3FDB5D1009E15CC0, 0x3C65B362CB974183), // 4.275550934302820943209668568888e-1
    F64x2::from_bits(0x3FD87DE2A6AEA963, 0xBC672CEDD3D5A610), // 3.826834323650897717284599840304e-1
    F64x2::from_bits(0x3FD58F9A75AB1FDD, 0xBC1EFDC0D58CF620), // 3.368898533922200506892532126191e-1
    F64x2::from_bits(0x3FD294062ED59F06, 0xBC75D28DA2C4612D), // 2.902846772544623676361923758174e-1
    F64x2::from_bits(0x3FCF19F97B215F1B, 0xBC642DEEF11DA2C4), // 2.429801799032638899482741620775e-1
    F64x2::from_bits(0x3FC8F8B83C69A60B, 0xBC626D19B9FF8D82), // 1.950903220161282678482848684770e-1
    F64x2::from_bits(0x3FC2C8106E8E613A, 0x3C513000A89A11E0), // 1.467304744553617516588501296467e-1
    F64x2::from_bits(0x3FB917A6BC29B42C, 0xBC3E2718D26ED688), // 9.801714032956060199419556388864e-2
    F64x2::from_bits(0x3FA91F65F10DD814, 0xBC2912BD0D569A90), // 4.906767432741801425495497694268e-2
    F64x2::from_bits(0x0000000000000000, 0x0000000000000000), // 0
    F64x2::from_bits(0xBFA91F65F10DD814, 0x3C2912BD0D569A90), // -4.906767432741801425495497694268e-2
    F64x2::from_bits(0xBFB917A6BC29B42C, 0x3C3E2718D26ED688), // -9.801714032956060199419556388864e-2
    F64x2::from_bits(0xBFC2C8106E8E613A, 0xBC513000A89A11E0), // -1.467304744553617516588501296467e-1
    F64x2::from_bits(0xBFC8F8B83C69A60B, 0x3C626D19B9FF8D82), // -1.950903220161282678482848684770e-1
    F64x2::from_bits(0xBFCF19F97B215F1B, 0x3C642DEEF11DA2C4), // -2.429801799032638899482741620775e-1
    F64x2::from_bits(0xBFD294062ED59F06, 0x3C75D28DA2C4612D), // -2.902846772544623676361923758174e-1
    F64x2::from_bits(0xBFD58F9A75AB1FDD, 0x3C1EFDC0D58CF620), // -3.368898533922200506892532126191e-1
    F64x2::from_bits(0xBFD87DE2A6AEA963, 0x3C672CEDD3D5A610), // -3.826834323650897717284599840304e-1
    F64x2::from_bits(0xBFDB5D1009E15CC0, 0xBC65B362CB974183), // -4.275550934302820943209668568888e-1
    F64x2::from_bits(0xBFDE2B5D3806F63B, 0xBC5E0D891D3C6841), // -4.713967368259976485563876259053e-1
    F64x2::from_bits(0xBFE073879922FFEE, 0x3C8A5A014347406C), // -5.141027441932217265936938389688e-1
    F64x2::from_bits(0xBFE1C73B39AE68C8, 0xBC8B25DD267F6600), // -5.555702330196022247428308139485e-1
    F64x2::from_bits(0xBFE30FF7FCE17035, 0x3C6EFCC626F74A6F), // -5.956993044924333434670365288300e-1
    F64x2::from_bits(0xBFE44CF325091DD6, 0xBC68076A2CFDC6B3), // -6.343932841636454982151716132255e-1
    F64x2::from_bits(0xBFE57D69348CECA0, 0x3C875720992BFBB2), // -6.715589548470184006253768504274e-1
    F64x2::from_bits(0xBFE6A09E667F3BCD, 0x3C8BDD3413B26456), // -7.071067811865475244008443621048e-1
    F64x2::from_bits(0xBFE7B5DF226AAFAF, 0x3C70F537ACDF0AD7), // -7.409511253549590911756168974952e-1
    F64x2::from_bits(0xBFE8BC806B151741, 0x3C82C5E12ED1336D), // -7.730104533627369608109066097585e-1
    F64x2::from_bits(0xBFE9B3E047F38741, 0x3C830EE286712474), // -8.032075314806449098066765129631e-1
    F64x2::from_bits(0xBFEA9B66290EA1A3, 0xBC39F630E8B6DAC8), // -8.314696123025452370787883776179e-1
    F64x2::from_bits(0xBFEB728345196E3E, 0x3C8BC69F324E6D61), // -8.577286100002720699022699842848e-1
    F64x2::from_bits(0xBFEC38B2F180BDB1, 0x3C76E0B1757C8D07), // -8.819212643483550297127568636604e-1
    F64x2::from_bits(0xBFECED7AF43CC773, 0x3C5E7B6BB5AB58AE), // -9.039892931234433315862002972305e-1
    F64x2::from_bits(0xBFED906BCF328D46, 0xBC7457E610231AC2), // -9.238795325112867561281831893968e-1
    F64x2::from_bits(0xBFEE212104F686E5, 0x3C8014C76C126527), // -9.415440651830207784125094025995e-1
    F64x2::from_bits(0xBFEE9F4156C62DDA, 0xBC8760B1E2E3F81E), // -9.569403357322088649357978869803e-1
    F64x2::from_bits(0xBFEF0A7EFB9230D7, 0xBC752C7ADC6B4989), // -9.700312531945439926039842072861e-1
    F64x2::from_bits(0xBFEF6297CFF75CB0, 0xBC7562172A361FD3), // -9.807852804032304491261822361342e-1
    F64x2::from_bits(0xBFEFA7557F08A517, 0x3C87A0A8CA13571F), // -9.891765099647809734516737380162e-1
    F64x2::from_bits(0xBFEFD88DA3D12526, 0x3C887DF6378811C7), // -9.951847266721968862448369531095e-1
    F64x2::from_bits(0xBFEFF621E3796D7E, 0x3C6C57BC2E24AA15), // -9.987954562051723927147716047591e-1
    F64x2::from_bits(0xBFF0000000000000, 0x0000000000000000), // -1.000000000000000000000000000000e0
    F64x2::from_bits(0xBFEFF621E3796D7E, 0x3C6C57BC2E24AA15), // -9.987954562051723927147716047591e-1
    F64x2::from_bits(0xBFEFD88DA3D12526, 0x3C887DF6378811C7), // -9.951847266721968862448369531095e-1
    F64x2::from_bits(0xBFEFA7557F08A517, 0x3C87A0A8CA13571F), // -9.891765099647809734516737380162e-1
    F64x2::from_bits(0xBFEF6297CFF75CB0, 0xBC7562172A361FD3), // -9.807852804032304491261822361342e-1
    F64x2::from_bits(0xBFEF0A7EFB9230D7, 0xBC752C7ADC6B4989), // -9.700312531945439926039842072861e-1
    F64x2::from_bits(0xBFEE9F4156C62DDA, 0xBC8760B1E2E3F81E), // -9.569403357322088649357978869803e-1
    F64x2::from_bits(0xBFEE212104F686E5, 0x3C8014C76C126527), // -9.415440651830207784125094025995e-1
    F64x2::from_bits(0xBFED906BCF328D46, 0xBC7457E610231AC2), // -9.238795325112867561281831893968e-1
    F64x2::from_bits(0xBFECED7AF43CC773, 0x3C5E7B6BB5AB58AE), // -9.039892931234433315862002972305e-1
    F64x2::from_bits(0xBFEC38B2F180BDB1, 0x3C76E0B1757C8D07), // -8.819212643483550297127568636604e-1
    F64x2::from_bits(0xBFEB728345196E3E, 0x3C8BC69F324E6D61), // -8.577286100002720699022699842848e-1
    F64x2::from_bits(0xBFEA9B66290EA1A3, 0xBC39F630E8B6DAC8), // -8.314696123025452370787883776179e-1
    F64x2::from_bits(0xBFE9B3E047F38741, 0x3C830EE286712474), // -8.032075314806449098066765129631e-1
    F64x2::from_bits(0xBFE8BC806B151741, 0x3C82C5E12ED1336D), // -7.730104533627369608109066097585e-1
    F64x2::from_bits(0xBFE7B5DF226AAFAF, 0x3C70F537ACDF0AD7), // -7.409511253549590911756168974952e-1
    F64x2::from_bits(0xBFE6A09E667F3BCD, 0x3C8BDD3413B26456), // -7.071067811865475244008443621048e-1
    F64x2::from_bits(0xBFE57D69348CECA0, 0x3C875720992BFBB2), // -6.715589548470184006253768504274e-1
    F64x2::from_bits(0xBFE44CF325091DD6, 0xBC68076A2CFDC6B3), // -6.343932841636454982151716132255e-1
    F64x2::from_bits(0xBFE30FF7FCE17035, 0x3C6EFCC626F74A6F), // -5.956993044924333434670365288300e-1
    F64x2::from_bits(0xBFE1C73B39AE68C8, 0xBC8B25DD267F6600), // -5.555702330196022247428308139485e-1
    F64x2::from_bits(0xBFE073879922FFEE, 0x3C8A5A014347406C), // -5.141027441932217265936938389688e-1
    F64x2::from_bits(0xBFDE2B5D3806F63B, 0xBC5E0D891D3C6841), // -4.713967368259976485563876259053e-1
    F64x2::from_bits(0xBFDB5D1009E15CC0, 0xBC65B362CB974183), // -4.275550934302820943209668568888e-1
    F64x2::from_bits(0xBFD87DE2A6AEA963, 0x3C672CEDD3D5A610), // -3.826834323650897717284599840304e-1
    F64x2::from_bits(0xBFD58F9A75AB1FDD, 0x3C1EFDC0D58CF620), // -3.368898533922200506892532126191e-1
    F64x2::from_bits(0xBFD294062ED59F06, 0x3C75D28DA2C4612D), // -2.902846772544623676361923758174e-1
    F64x2::from_bits(0xBFCF19F97B215F1B, 0x3C642DEEF11DA2C4), // -2.429801799032638899482741620775e-1
    F64x2::from_bits(0xBFC8F8B83C69A60B, 0x3C626D19B9FF8D82), // -1.950903220161282678482848684770e-1
    F64x2::from_bits(0xBFC2C8106E8E613A, 0xBC513000A89A11E0), // -1.467304744553617516588501296467e-1
    F64x2::from_bits(0xBFB917A6BC29B42C, 0x3C3E2718D26ED688), // -9.801714032956060199419556388864e-2
    F64x2::from_bits(0xBFA91F65F10DD814, 0x3C2912BD0D569A90), // -4.906767432741801425495497694268e-2
];

impl Reduced {
    /// Returns `sin(n * π/64 + b)` as a non-normalized sum `hi + lo`.
    #[inline]
    fn sin_k(&self, n: u8) -> Sum2 {
        let s = SIN_PI_TBL[usize::from(n & 127)];
        let c = SIN_PI_TBL[usize::from(n.wrapping_add(32) & 127)];
        self.sin_a(s, c)
    }

    /// Returns `cos(n * π/64 + b)` as a non-normalized sum `hi + lo`.
    #[inline]
    fn cos_k(&self, n: u8) -> Sum2 {
        self.sin_k(n.wrapping_add(32))
    }

    /// Returns `tan(n * π/64 + b)`.
    #[inline]
    fn tan_k(&self, n: u8) -> f64 {
        let s = self.sin_k(n).normalize();
        let c = self.cos_k(n).normalize();
        if c.hi() == 0.0 {
            // Only when `b` is zero and `n` is an odd multiple of 32, that
            // is, an exact odd multiple of a right angle in degrees or half
            // revolutions: the result is infinite.
            return s.hi() / c.hi();
        }

        // q0 ~= s / c
        let inv = 1.0 / c.hi();
        let q0 = (s.hi() * inv).purify();
        // Exact product
        let p = F64x2::mul11(q0, c.hi());
        // (s - q0 * c), where `s.hi() - p.hi()` is exact by Sterbenz lemma
        let rem = ((s.hi() - p.hi()) - p.lo()) + (s.lo() - q0 * c.lo());
        let q1 = rem * inv;

        (q0 + q1).purify()
    }

    /// Returns `sin(a + b)` as a non-normalized sum `hi + lo`, where
    /// `s = sin(a)` and `c = cos(a)`.
    ///
    /// Requires `|s.hi()| >= |c.hi() * b|` or `s.hi() == 0`.
    #[inline(always)]
    fn sin_a(&self, s: F64x2, c: F64x2) -> Sum2 {
        // GENERATE: sin_poly f64 3 0.0247
        const K3: f64 = f64::from_bits(0xBFC5555555555551); // -1.6666666666666655e-1
        const K5: f64 = f64::from_bits(0x3F8111111107054E); // 8.333333332191242e-3
        const K7: f64 = f64::from_bits(0xBF2A01845C6DE673); // -1.984094689893471e-4

        // GENERATE: cos_poly f64 3 0.0247
        const K4: f64 = f64::from_bits(0x3FA5555555555553); // 4.166666666666665e-2
        const K6: f64 = f64::from_bits(0xBF56C16C16B6D991); // -1.3888888887386467e-3
        const K8: f64 = f64::from_bits(0x3EFA01873EEED8D9); // 2.4801225602743462e-5

        let (sh, sl) = (s.hi(), s.lo());
        let (ch, cl) = (c.hi(), c.lo());
        let bh = self.bh;
        let bl = self.bl;

        // With b = bh + bl, neglecting O(bl^2) terms:
        //   sin(a + b) = sin(a + bh) + bl * cos(a + bh)
        //   sin(a + bh) = s + c * bh + s * (cos(bh) - 1) + c * (sin(bh) - bh)
        //   cos(a + bh) = c * (1 - bh^2 / 2) - s * bh + O(bh^3)
        let b2 = bh * bh;
        let hb2 = 0.5 * b2;
        // cos(bh) - 1 = -hb2 + cm1_rem
        let cm1_rem = b2 * horner!(b2, b2, [K4, K6, K8]);
        // sin(bh) - bh
        let smb = bh * horner!(b2, b2, [K3, K5, K7]);

        // c * bh (exact)
        let p = F64x2::mul11(ch, bh);
        // s + c * bh (exact)
        let r = F64x2::fast_add11(sh, p.hi());

        // The remaining terms, from smallest to largest. The largest one,
        // `(sh + sl) * hb2`, is added last, so its error is only affected by
        // three roundings (`b2`, `sh * hb2` and `lo`).
        let t = (r.lo() + (sl + p.lo()))
            + (cl * bh + bl * (ch - (sh * bh + ch * hb2)))
            + ((sh * cm1_rem - sl * hb2) + ch * smb);
        let lo = t - sh * hb2;

        Sum2 { hi: r.hi(), lo }
    }
}

/// Non-normalized sum `hi + lo`, with `|lo|` much smaller than `|hi|`.
struct Sum2 {
    hi: f64,
    lo: f64,
}

impl Sum2 {
    #[inline]
    fn to_f64(&self) -> f64 {
        (self.hi + self.lo).purify()
    }

    #[inline]
    fn normalize(&self) -> F64x2 {
        F64x2::fast_add11(self.hi, self.lo.purify())
    }
}
