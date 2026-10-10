mod cbrt;
mod exp;
mod f64x2;
mod gamma;
mod hyperbolic;
mod hypot;
mod inv_hyperbolic;
mod inv_trigonometric;
mod log;
mod pow;
mod sqrt;
mod trigonometric;

pub(crate) use exp::{EXP2_TBL_BITS, exp2_tbl};
pub(crate) use gamma::ln_gamma_near_zero;
pub(crate) use inv_hyperbolic::{acosh_large_corr, asinh_large_corr};
pub(crate) use inv_trigonometric::{atan_index, atan_tbl};
pub(crate) use log::{ln_1p_q, ln_tbl, split_ln_arg};
pub(crate) use trigonometric::{FRAC_64_PI, reduce_rad_large_sum, sin_cos_pi_128};

use crate::traits::Float as _;

impl crate::traits::Float for f64 {
    type Raw = u64;

    type RawExp = u16;

    type Exp = i16;

    const BITS: u8 = 64;
    const MANT_BITS: u8 = 52;
    const EXP_BITS: u8 = 11;

    const SIGN_MASK: Self::Raw = 1 << (<Self as crate::traits::Float>::BITS - 1);
    const EXP_MASK: Self::Raw = ((1 << Self::EXP_BITS) - 1) << Self::MANT_BITS;
    const MANT_MASK: Self::Raw = (1 << Self::MANT_BITS) - 1;

    const EXP_OFFSET: Self::RawExp = (1 << (Self::EXP_BITS - 1)) - 1;
    const MAX_RAW_EXP: Self::RawExp = (Self::EXP_MASK >> Self::MANT_BITS) as Self::RawExp;

    const MIN_NORMAL_EXP: Self::Exp = -<Self as crate::traits::Float>::MAX_EXP + 1;
    const MAX_EXP: Self::Exp = (Self::MAX_RAW_EXP >> 1) as Self::Exp;

    const INFINITY: Self = Self::INFINITY;
    const NEG_INFINITY: Self = Self::NEG_INFINITY;
    const NAN: Self = Self::NAN;

    const ZERO: Self = 0.0;
    const HALF: Self = 0.5;
    const ONE: Self = 1.0;
    const TWO: Self = 2.0;

    #[cfg(test)]
    const LARGEST: Self = Self::MAX;

    #[inline]
    fn purify(self) -> Self {
        if cfg!(all(
            any(target_arch = "x86", target_arch = "x86_64"),
            not(target_feature = "sse2")
        )) {
            // Workaround X87 rounding issues
            // `read_volatile` documentation says "Volatile operations are intended
            // to act on I/O memory, and are guaranteed to not be elided or...". Not
            // being elided means not being optimized away. Using `read_volatile::<f64>`
            // guarantees that the returned value is the result of a 8-byte memory read,
            // so it cannot have precision beyond a `f64`.
            unsafe { core::ptr::read_volatile(&self) }
        } else {
            self
        }
    }

    #[inline]
    fn to_raw(self) -> Self::Raw {
        self.to_bits()
    }

    #[inline]
    fn from_raw(raw: Self::Raw) -> Self {
        Self::from_bits(raw)
    }

    #[inline]
    fn raw_exp_to_exp(e: Self::RawExp) -> Self::Exp {
        e.wrapping_sub(Self::EXP_OFFSET) as i16
    }

    #[inline]
    fn exp_to_raw_exp(e: Self::Exp) -> Self::RawExp {
        (e as Self::RawExp).wrapping_add(Self::EXP_OFFSET)
    }

    #[inline]
    fn abs(self) -> Self {
        self.abs()
    }

    #[inline]
    fn copysign(self, y: Self) -> Self {
        self.copysign(y)
    }

    #[cfg(test)]
    fn parse(s: &str) -> Self {
        s.parse().unwrap()
    }
}

impl crate::sealed::SealedMath for f64 {}

impl crate::FloatMath for f64 {
    fn abs(x: Self) -> Self {
        x.abs()
    }

    fn copysign(x: Self, y: Self) -> Self {
        x.copysign(y)
    }

    fn round(x: Self) -> Self {
        crate::generic::round(x)
    }

    fn trunc(x: Self) -> Self {
        crate::generic::trunc(x)
    }

    fn ceil(x: Self) -> Self {
        crate::generic::ceil(x)
    }

    fn floor(x: Self) -> Self {
        crate::generic::floor(x)
    }

    fn scalbn(x: Self, y: i32) -> Self {
        crate::generic::scalbn(x, y)
    }

    fn frexp(x: Self) -> (Self, i32) {
        crate::generic::frexp(x)
    }

    fn hypot(x: Self, y: Self) -> Self {
        crate::generic::hypot(x, y)
    }

    fn sqrt(x: Self) -> Self {
        crate::generic::sqrt(x)
    }

    fn cbrt(x: Self) -> Self {
        crate::generic::cbrt(x)
    }

    fn exp(x: Self) -> Self {
        crate::generic::exp(x)
    }

    fn exp_m1(x: Self) -> Self {
        crate::generic::exp_m1(x)
    }

    fn exp2(x: Self) -> Self {
        crate::generic::exp2(x)
    }

    fn exp10(x: Self) -> Self {
        crate::generic::exp10(x)
    }

    fn ln(x: Self) -> Self {
        crate::generic::ln(x)
    }

    fn ln_1p(x: Self) -> Self {
        crate::generic::ln_1p(x)
    }

    fn log2(x: Self) -> Self {
        crate::generic::log2(x)
    }

    fn log10(x: Self) -> Self {
        crate::generic::log10(x)
    }

    fn pow(x: Self, y: Self) -> Self {
        crate::generic::pow(x, y)
    }

    fn powi(x: Self, y: i32) -> Self {
        crate::generic::powi(x, y)
    }

    fn sin(x: Self) -> Self {
        crate::generic::sin(x)
    }

    fn cos(x: Self) -> Self {
        crate::generic::cos(x)
    }

    fn sin_cos(x: Self) -> (Self, Self) {
        crate::generic::sin_cos(x)
    }

    fn tan(x: Self) -> Self {
        crate::generic::tan(x)
    }

    fn sind(x: Self) -> Self {
        crate::generic::sind(x)
    }

    fn cosd(x: Self) -> Self {
        crate::generic::cosd(x)
    }

    fn sind_cosd(x: Self) -> (Self, Self) {
        crate::generic::sind_cosd(x)
    }

    fn tand(x: Self) -> Self {
        crate::generic::tand(x)
    }

    fn sinpi(x: Self) -> Self {
        crate::generic::sinpi(x)
    }

    fn cospi(x: Self) -> Self {
        crate::generic::cospi(x)
    }

    fn sinpi_cospi(x: Self) -> (Self, Self) {
        crate::generic::sinpi_cospi(x)
    }

    fn tanpi(x: Self) -> Self {
        crate::generic::tanpi(x)
    }

    fn asin(x: Self) -> Self {
        crate::generic::asin(x)
    }

    fn acos(x: Self) -> Self {
        crate::generic::acos(x)
    }

    fn atan(x: Self) -> Self {
        crate::generic::atan(x)
    }

    fn atan2(y: Self, x: Self) -> Self {
        crate::generic::atan2(y, x)
    }

    fn asind(x: Self) -> Self {
        crate::generic::asind(x)
    }

    fn acosd(x: Self) -> Self {
        crate::generic::acosd(x)
    }

    fn atand(x: Self) -> Self {
        crate::generic::atand(x)
    }

    fn atan2d(y: Self, x: Self) -> Self {
        crate::generic::atan2d(y, x)
    }

    fn asinpi(x: Self) -> Self {
        crate::generic::asinpi(x)
    }

    fn acospi(x: Self) -> Self {
        crate::generic::acospi(x)
    }

    fn atanpi(x: Self) -> Self {
        crate::generic::atanpi(x)
    }

    fn atan2pi(y: Self, x: Self) -> Self {
        crate::generic::atan2pi(y, x)
    }

    fn sinh(x: Self) -> Self {
        crate::generic::sinh(x)
    }

    fn cosh(x: Self) -> Self {
        crate::generic::cosh(x)
    }

    fn sinh_cosh(x: Self) -> (Self, Self) {
        crate::generic::sinh_cosh(x)
    }

    fn tanh(x: Self) -> Self {
        crate::generic::tanh(x)
    }

    fn asinh(x: Self) -> Self {
        crate::generic::asinh(x)
    }

    fn acosh(x: Self) -> Self {
        crate::generic::acosh(x)
    }

    fn atanh(x: Self) -> Self {
        crate::generic::atanh(x)
    }

    fn gamma(x: Self) -> Self {
        crate::generic::gamma(x)
    }

    fn ln_gamma(x: Self) -> (Self, i8) {
        crate::generic::ln_gamma(x)
    }
}

/// Returns approximations of `(sqrt(x), 1 / sqrt(x))`, with relative errors
/// of at most about 2^-50.7 and 2^-51, for `2^-1021 <= x < 2^1022`.
///
/// Outside that range, `x / 2` or `r^2` (see below) is subnormal, which
/// increases the errors (up to about 2^-50.3) and is slow.
#[inline]
pub(crate) fn fast_sqrt(x: f64) -> (f64, f64) {
    // Initial approximation of r ~= 1 / sqrt(x) using bit fiddling, with a
    // relative error less than 2^-4.8.
    let hx = 0.5 * x;
    let mut r = f64::from_bits(0x5FE6_EB50_C7B5_37A9_u64.wrapping_sub(x.to_bits() >> 1));
    // Four Newton iterations, each one roughly doubling the number of correct
    // bits:
    // r' = r * (3 - x * r^2) / 2 = 1.5 * r - (x / 2 * r) * r^2
    // where the last form has a chain of three dependent operations instead
    // of four, but one more rounding, so the last iteration leaves an error
    // of up to 4 * 2^-53 instead of 3 * 2^-53 (and the final product adds
    // 2^-53 to the error of the square root).
    for _ in 0..4 {
        r = 1.5 * r - (hx * r) * (r * r);
    }
    (r * x, r)
}

/// Returns `(p, e)` such that `p` is `x^2` rounded to `f64` and `p + e ~= x^2`,
/// with an error less than 2^-104 * x^2 when `x^2 >= 2^-960` (so a subnormal
/// `b^2` does not lose too much accuracy), and no error when
/// `1 <= x < 1 + 2^-26` (as `inv_hyperbolic::acosh_parts` requires).
#[inline]
fn square_parts(x: f64) -> (f64, f64) {
    // x = a + b, where `a` has 26 bits and `b` at most 27, so `a^2` and
    // `2 * a * b` are exact, and so is `a^2 - p` (by Sterbenz lemma). Only
    // `b^2` (less than 2^-50 * x^2) and the following sums can be rounded.
    // When `1 <= x < 1 + 2^-26`, `a = 1` and `b` has at most 26 bits, so
    // `b^2` is exact, and so are the sums, whose results are multiples of
    // 2^-104 with magnitudes below 2^-50.
    let a = x.split_hi();
    let b = x - a;
    let p = (x * x).purify();
    let e = ((a * a - p) + 2.0 * a * b) + b * b;
    (p, e)
}

/// Returns `(hi, lo)` such that `hi + lo ~= sqrt(w_hi + w_lo)`, with a
/// relative error of about 2^-100, where `2^-960 <= w_hi < 2^1022` and
/// `|w_lo| <= 2^-51 * w_hi`.
#[inline]
fn sqrt_parts(w_hi: f64, w_lo: f64) -> (f64, f64) {
    // y ~= sqrt(w_hi) and r ~= 1 / sqrt(w_hi), with relative errors of at
    // most about 2^-50.7 and 2^-51 (see `fast_sqrt`), which the Newton
    // iteration below roughly squares
    let (y, r) = fast_sqrt(w_hi);
    let y = y.purify();
    // One Newton iteration, with the residual calculated accurately:
    // sqrt(w) ~= y + (w - y^2) / (2 * y) ~= y + (w - y^2) * r / 2
    // `w_hi - p` is exact (by Sterbenz lemma).
    let (p, pe) = square_parts(y);
    let res = ((w_hi - p) - pe) + w_lo;
    (y, res * (0.5 * r))
}

/// Rounds `x` to the nearest integer, returning it as `f64` and its lowest
/// 8 bits (in two's complement).
///
/// `|x|` must be less than 2^51.
#[inline]
pub(crate) fn round_u8(x: f64) -> (f64, u8) {
    let (r, bits) = round_to_bits(x);
    (r, bits as u8)
}

/// Rounds `x` to the nearest integer, returning it as `f64` and `i32`.
///
/// `|x|` must be less than 2^31.
#[inline]
pub(crate) fn round_i32(x: f64) -> (f64, i32) {
    let (r, bits) = round_to_bits(x);
    (r, bits as i32)
}

/// Rounds `x` to the nearest integer, returning it as `f64` and a `u64`
/// whose lowest bits are those of the integer (in two's complement).
///
/// `|x|` must be less than 2^51.
#[inline]
fn round_to_bits(x: f64) -> (f64, u64) {
    // 1.5 * 2^52, adding it rounds to an integer, keeping the integer in the
    // lowest bits of the mantissa.
    const MAGIC: f64 = (3u64 << 51) as f64;

    let t = (x + MAGIC).purify();
    ((t - MAGIC).purify(), t.to_bits())
}

#[cfg(test)]
mod tests {
    use super::round_i32;
    use crate::traits::Float as _;

    #[test]
    fn test_exp2i_fast() {
        for e in -1022..=1023 {
            let x = f64::exp2i_fast(e);
            assert_eq!(x, f64::exp2(e as f64));
            assert_eq!(x.exponent(), e);
            assert_eq!(x.to_bits() & f64::MANT_MASK, 0);
        }
    }

    #[test]
    fn test_normalize_arg() {
        fn check(x: f64) {
            let (y, edelta) = x.normalize_arg();
            if x.raw_exp() == 0 {
                assert_eq!(y.to_bits(), (x * f64::exp2i_fast(52)).to_bits());
                assert_eq!(edelta, -52);
            } else {
                assert_eq!(y.to_bits(), x.to_bits());
                assert_eq!(edelta, 0);
            }
        }

        for sign in [0, 1 << 63] {
            check(f64::from_bits(sign));
            for i in 0..52 {
                check(f64::from_bits(sign | (1 << i)));
                check(f64::from_bits(sign | ((1 << (i + 1)) - 1)));
            }
            check(f64::from_bits(sign | 0x000F_EDCB_A987_6543));
            check(f64::from_bits(sign | f64::MIN_POSITIVE.to_bits()));
            check(f64::from_bits(sign | 1.5f64.to_bits()));
        }
    }

    #[test]
    fn test_round_i32() {
        for i in -300_000..=300_000 {
            let x = (f64::from(i) * 0.37).purify();
            let (xf, xi) = round_i32(x);
            assert_eq!(f64::from(xi), xf);
            // On x87, double rounding can move the result to the other
            // neighbouring integer.
            assert!((x - xf).abs() <= 0.5 + 1.0 / 2048.0, "x = {x:e}");

            // ties to even
            let x = f64::from(i) + 0.5;
            let (xf, xi) = round_i32(x);
            assert_eq!(f64::from(xi), xf);
            assert_eq!(xi, i + (i & 1));
        }
    }
}
