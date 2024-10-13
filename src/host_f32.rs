use crate::traits;

impl traits::Float for f32 {
    type Like = crate::f32::LikeF32;

    type Raw = u32;

    type RawExp = u8;

    type Exp = i16;

    const BITS: u8 = 32;
    const MANT_BITS: u8 = 23;
    const EXP_BITS: u8 = 8;

    const SIGN_MASK: Self::Raw = 1 << (Self::BITS - 1);
    const EXP_MASK: Self::Raw = ((1 << Self::EXP_BITS) - 1) << Self::MANT_BITS;
    const MANT_MASK: Self::Raw = (1 << Self::MANT_BITS) - 1;

    const EXP_OFFSET: Self::RawExp = (1 << (Self::EXP_BITS - 1)) - 1;
    const MAX_RAW_EXP: Self::RawExp = (Self::EXP_MASK >> Self::MANT_BITS) as Self::RawExp;

    const MIN_NORMAL_EXP: Self::Exp = -<Self as traits::Float>::MAX_EXP + 1;
    const MAX_EXP: Self::Exp = (Self::MAX_RAW_EXP >> 1) as Self::Exp;

    const INFINITY: Self = Self::INFINITY;

    #[inline]
    fn neg_infinity() -> Self {
        Self::NEG_INFINITY
    }

    const NAN: Self = Self::NAN;

    const ZERO: Self = 0.0;

    #[inline]
    fn half() -> Self {
        0.5
    }

    #[inline]
    fn one() -> Self {
        1.0
    }

    #[inline]
    fn two() -> Self {
        2.0
    }

    #[cfg(test)]
    #[inline]
    fn largest() -> Self {
        Self::MAX
    }

    #[inline]
    fn purify(self) -> Self {
        if cfg!(all(
            any(target_arch = "x86", target_arch = "x86_64"),
            not(target_feature = "sse2")
        )) {
            // Workaround X87 rounding issues
            // `read_volatile` documentation says "Volatile operations are intended
            // to act on I/O memory, and are guaranteed to not be elided or...". Not
            // being elided means not being optimized away. Using `read_volatile::<f32>`
            // guarantees that the returned value is the result of a 4-byte memory read,
            // so it cannot have precision beyond a `f32`.
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
        i16::from(e.wrapping_sub(Self::EXP_OFFSET) as i8)
    }

    #[inline]
    fn exp_to_raw_exp(e: Self::Exp) -> Self::RawExp {
        (e as Self::RawExp).wrapping_add(Self::EXP_OFFSET)
    }

    #[cfg(test)]
    #[inline]
    fn is_nan(self) -> bool {
        self.is_nan()
    }

    #[cfg(test)]
    fn parse(s: &str) -> Self {
        s.parse().unwrap()
    }
}

impl crate::f32::F32Like for f32 {}

impl crate::sealed::SealedMath for f32 {}

impl crate::FloatMath for f32 {
    fn abs(x: Self) -> Self {
        traits::Float::abs(x)
    }

    fn copysign(x: Self, y: Self) -> Self {
        traits::Float::copysign(x, y)
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

    fn log(x: Self) -> Self {
        crate::generic::log(x)
    }

    fn log_1p(x: Self) -> Self {
        crate::generic::log_1p(x)
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

    fn tgamma(x: Self) -> Self {
        crate::generic::tgamma(x)
    }

    fn lgamma(x: Self) -> (Self, i8) {
        crate::generic::lgamma(x)
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::Float as _;

    #[test]
    fn test_exp2i_fast() {
        for e in -126..=127 {
            let x = f32::exp2i_fast(e);
            assert_eq!(x, f32::exp2(e as f32));
            assert_eq!(x.exponent(), e);
            assert_eq!(x.to_bits() & f32::MANT_MASK, 0);
        }
    }
}
