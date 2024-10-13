use rustc_apfloat::Float as _;

use crate::traits;

/// Software single-precision floating-point number.
///
/// Only available when the `soft-float` feature is enabled.
#[derive(Copy, Clone, Default, PartialEq, PartialOrd)]
pub struct SoftF32(rustc_apfloat::ieee::Single);

impl SoftF32 {
    /// Returns infinity (∞).
    #[inline]
    pub fn infinity() -> Self {
        Self(rustc_apfloat::Float::INFINITY)
    }

    /// Returns negative infinity (-∞).
    #[inline]
    pub fn neg_infinity() -> Self {
        -Self(rustc_apfloat::Float::INFINITY)
    }

    /// Returns a Not a Number (NaN).
    #[inline]
    pub fn nan() -> Self {
        Self(rustc_apfloat::Float::NAN)
    }

    /// Returns the raw representation of the floating-point number.
    #[inline]
    pub fn to_bits(self) -> u32 {
        self.0.to_bits() as _
    }

    /// Creates a floating-point number from its raw representation.
    #[inline]
    pub fn from_bits(bits: u32) -> Self {
        Self(rustc_apfloat::Float::from_bits(bits.into()))
    }

    /// Converts the floating-point number to the native type.
    #[inline]
    pub fn to_host(self) -> f32 {
        f32::from_bits(self.to_bits())
    }

    /// Creates a soft-float from the native type.
    #[inline]
    pub fn from_host(value: f32) -> Self {
        Self::from_bits(value.to_bits())
    }

    /// Returns the floating point category of the number.
    #[inline]
    pub fn classify(self) -> core::num::FpCategory {
        match self.0.category() {
            rustc_apfloat::Category::Infinity => core::num::FpCategory::Infinite,
            rustc_apfloat::Category::NaN => core::num::FpCategory::Nan,
            rustc_apfloat::Category::Normal => {
                if self.0.is_denormal() {
                    core::num::FpCategory::Subnormal
                } else {
                    core::num::FpCategory::Normal
                }
            }
            rustc_apfloat::Category::Zero => core::num::FpCategory::Zero,
        }
    }

    /// Returns `true`` if the value is neither infinite nor NaN.
    #[inline]
    pub fn is_finite(self) -> bool {
        self.0.is_finite()
    }

    /// Returns `true` if the value is infinite.
    #[inline]
    pub fn is_infinite(self) -> bool {
        self.0.is_infinite()
    }

    /// Returns `true` if the value is NaN.
    #[inline]
    pub fn is_nan(self) -> bool {
        self.0.is_nan()
    }

    #[inline]
    fn round_to_integral(self, round: rustc_apfloat::Round) -> Self {
        Self(self.0.round_to_integral(round).value)
    }

    #[inline]
    fn scalbn(self, y: i32) -> Self {
        Self(self.0.scalbn(y))
    }
}

impl core::ops::Neg for SoftF32 {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self(-self.0)
    }
}

impl core::ops::Add for SoftF32 {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self((self.0 + rhs.0).value)
    }
}

impl core::ops::Sub for SoftF32 {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self((self.0 - rhs.0).value)
    }
}

impl core::ops::Mul for SoftF32 {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self((self.0 * rhs.0).value)
    }
}

impl core::ops::Div for SoftF32 {
    type Output = Self;

    #[inline]
    fn div(self, rhs: Self) -> Self {
        Self((self.0 / rhs.0).value)
    }
}

impl core::fmt::Debug for SoftF32 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        core::fmt::Debug::fmt(&self.0, f)
    }
}

impl core::fmt::Display for SoftF32 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        core::fmt::Display::fmt(&self.0, f)
    }
}

impl traits::CastFrom<u8> for SoftF32 {
    #[inline]
    fn cast_from(value: u8) -> Self {
        Self(rustc_apfloat::Float::from_i128(value.into()).value)
    }
}

impl traits::CastFrom<i32> for SoftF32 {
    #[inline]
    fn cast_from(value: i32) -> Self {
        Self(rustc_apfloat::Float::from_i128(value.into()).value)
    }
}

impl traits::CastFrom<u32> for SoftF32 {
    #[inline]
    fn cast_from(value: u32) -> Self {
        Self(rustc_apfloat::Float::from_u128(value.into()).value)
    }
}

impl traits::CastFrom<i64> for SoftF32 {
    #[inline]
    fn cast_from(value: i64) -> Self {
        Self(rustc_apfloat::Float::from_i128(value.into()).value)
    }
}

impl traits::CastFrom<u64> for SoftF32 {
    #[inline]
    fn cast_from(value: u64) -> Self {
        Self(rustc_apfloat::Float::from_u128(value.into()).value)
    }
}

impl traits::CastFrom<i16> for SoftF32 {
    #[inline]
    fn cast_from(value: i16) -> Self {
        Self(rustc_apfloat::Float::from_i128(value.into()).value)
    }
}

impl traits::Float for SoftF32 {
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

    const INFINITY: Self = Self(rustc_apfloat::Float::INFINITY);

    #[inline]
    fn neg_infinity() -> Self {
        -Self(rustc_apfloat::Float::INFINITY)
    }

    const NAN: Self = Self(rustc_apfloat::Float::NAN);

    const ZERO: Self = Self(rustc_apfloat::Float::ZERO);

    #[inline]
    fn half() -> Self {
        Self::from_raw(0x3F000000)
    }

    #[inline]
    fn one() -> Self {
        Self::from_raw(0x3F800000)
    }

    #[inline]
    fn two() -> Self {
        Self::from_raw(0x40000000)
    }

    #[cfg(test)]
    #[inline]
    fn largest() -> Self {
        Self(rustc_apfloat::Float::largest())
    }

    #[inline]
    fn purify(self) -> Self {
        self
    }

    #[inline]
    fn to_raw(self) -> Self::Raw {
        self.0.to_bits() as _
    }

    #[inline]
    fn from_raw(raw: Self::Raw) -> Self {
        Self(rustc_apfloat::Float::from_bits(raw.into()))
    }

    #[inline]
    fn raw_exp_to_exp(e: Self::RawExp) -> Self::Exp {
        i16::from(e.wrapping_sub(Self::EXP_OFFSET) as i8)
    }

    #[inline]
    fn exp_to_raw_exp(e: Self::Exp) -> Self::RawExp {
        (e as Self::RawExp).wrapping_add(Self::EXP_OFFSET)
    }

    #[inline]
    fn sign(self) -> bool {
        self.0.is_negative()
    }

    #[cfg(test)]
    #[inline]
    fn is_nan(self) -> bool {
        self.0.is_nan()
    }

    #[inline]
    fn abs(self) -> Self {
        Self(self.0.abs())
    }

    #[inline]
    fn copysign(self, y: Self) -> Self {
        Self(self.0.copy_sign(y.0))
    }

    #[cfg(test)]
    fn parse(s: &str) -> Self {
        Self(s.parse().unwrap())
    }
}

impl crate::f32::F32Like for SoftF32 {}

impl crate::sealed::SealedMath for SoftF32 {}

impl crate::FloatMath for SoftF32 {
    fn abs(x: Self) -> Self {
        Self(x.0.abs())
    }

    fn copysign(x: Self, y: Self) -> Self {
        Self(x.0.copy_sign(y.0))
    }

    fn round(x: Self) -> Self {
        x.round_to_integral(rustc_apfloat::Round::NearestTiesToAway)
    }

    fn trunc(x: Self) -> Self {
        x.round_to_integral(rustc_apfloat::Round::TowardZero)
    }

    fn ceil(x: Self) -> Self {
        x.round_to_integral(rustc_apfloat::Round::TowardPositive)
    }

    fn floor(x: Self) -> Self {
        x.round_to_integral(rustc_apfloat::Round::TowardNegative)
    }

    fn scalbn(x: Self, y: i32) -> Self {
        x.scalbn(y)
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
