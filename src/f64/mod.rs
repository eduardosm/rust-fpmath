mod cbrt;
mod exp;
mod f64x2;
mod gamma;
mod hyperbolic;
mod hypot;
mod inv_hyperbolic;
mod inv_trigonometric;
mod log;
mod log_core;
mod pow;
mod sqrt;
mod trigonometric;

pub(crate) use exp::{EXP2_TBL_BITS, exp2_tbl};
pub(crate) use inv_hyperbolic::{acosh_large_corr, asinh_large_corr};
pub(crate) use log::{ln_tbl, split_ln_arg};
pub(crate) use trigonometric::{FRAC_64_PI, reduce_rad_large_sum, sin_cos_pi_64};

use crate::traits::Float as _;

impl crate::traits::Float for f64 {
    type Raw = u64;

    type SRaw = i64;

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

/// Return an approximation of `(sqrt(x), 1/sqrt(x))`
pub(crate) fn fast_sqrt(x: f64) -> (f64, f64) {
    // Initial approximation of r ~= 1/sqrt(x) using bit fiddling.
    let hx = 0.5 * x;
    let mut r = f64::from_bits(0x5FE6_EB50_C7B5_37A9_u64.wrapping_sub(x.to_bits() >> 1));
    // Four Newton iterations to reduce error to nearly 1 ULP
    // r_next = r * (3 - hi * r^2) / 2
    for _ in 0..4 {
        r *= 1.5 - hx * r * r;
    }
    (r * x, r)
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
