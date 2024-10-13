#![warn(
    rust_2018_idioms,
    trivial_casts,
    trivial_numeric_casts,
    unreachable_pub,
    unused_qualifications
)]
#![no_std]

//! A pure-Rust floating point math library with support for native floating
//! point and soft-floats.
//!
//! This crate is `no_std`.
//!
//! The following math functions are implemented:
//!
//! * Sign operations ([`abs`], [`copysign`]).
//! * Rounding ([`round`], [`trunc`], [`ceil`], [`floor`]).
//! * Exponential ([`exp`], [`exp_m1`], [`exp2`], [`exp10`]).
//! * Logarithmic ([`log`], [`log_1p`], [`log2`], [`log10`]).
//! * Power ([`pow`], [`powi`]).
//! * Trigonometric
//!   - Radians ([`sin`], [`cos`], [`sin_cos`], [`tan`]).
//!   - Degrees ([`sind`], [`cosd`], [`sind_cosd`], [`tand`]).
//!   - Half-revolutions ([`sinpi`], [`cospi`], [`sinpi_cospi`], [`tanpi`]).
//! * Inverse trigonometric
//!   - Radians ([`asin`], [`acos`], [`atan`], [`atan2`]).
//!   - Degrees ([`asind`], [`acosd`], [`atand`], [`atan2d`]).
//!   - Half-revolutions ([`asinpi`], [`acospi`], [`atanpi`], [`atan2pi`]).
//! * Hyperbolic ([`sinh`], [`cosh`], [`sinh_cosh`], [`tanh`]).
//! * Inverse hyperbolic ([`asinh`], [`acosh`], [`atanh`]).
//! * Gamma ([`tgamma`], [`lgamma`]).
//!
//! All functions are implemted for the native floating point types [`prim@f32`]
//! and [`prim@f64`].
//!
//! The soft-float types [`SoftF32`] and [`SoftF64`] are also provided. They
//! also support all the above functions and provide consistent bit-to-bit
//! behavior across platforms. They are available when the `soft-float` feature
//! is enabled (disabled by default).
//!
//! The [`FloatMath`] trait is used to identify types that support the math
//! functions.

// TODO:
// * Error function and complementary (erf, erfc)
// * Bessel functions (j0, y0, j1, y1, jn, yn)

// Uncomment to use `dbg!`
//extern crate std;

macro_rules! horner {
    ($outer_x:ident, $inner_x:ident, [$coef:expr]) => {
        $outer_x * $coef
    };
    ($outer_x:ident, $inner_x:ident, [$coef0:expr, $($coefs:expr),+]) => {
        $outer_x * ($coef0 + horner!($inner_x, $inner_x, [$($coefs),+]))
    };
}

#[cfg(test)]
macro_rules! assert_is_nan {
    ($value:expr) => {{
        let value = $value;
        if !$crate::traits::Float::is_nan(value) {
            panic!("assertion failed: `({value:?}).is_nan()`");
        }
    }};
}

#[cfg(test)]
macro_rules! assert_total_eq {
    ($lhs:expr, $rhs:expr) => {
        let [lhs, rhs] = [$lhs, $rhs];
        let lhs = $crate::traits::Float::purify(lhs);
        let rhs = $crate::traits::Float::purify(rhs);
        let lhs_raw = $crate::traits::Float::to_raw(lhs);
        let rhs_raw = $crate::traits::Float::to_raw(rhs);
        if lhs_raw != rhs_raw {
            panic!("assertion failed: `{lhs:?} == {rhs:?}` (using totalOrder)");
        }
    };
}

mod double;
mod f32;
mod f64;
mod generic;
mod host_f32;
mod host_f64;
mod int;
#[cfg(feature = "soft-float")]
mod soft_f32;
#[cfg(feature = "soft-float")]
mod soft_f64;
mod traits;

#[cfg(feature = "soft-float")]
pub use soft_f32::SoftF32;
#[cfg(feature = "soft-float")]
pub use soft_f64::SoftF64;

mod sealed {
    pub trait SealedMath {}
}

/// Floating point types with math functions.
pub trait FloatMath: sealed::SealedMath + Sized {
    /// See the [`abs`] function.
    fn abs(x: Self) -> Self;

    /// See the [`copysign`] function.
    fn copysign(x: Self, y: Self) -> Self;

    /// See the [`round`] function.
    fn round(x: Self) -> Self;

    /// See the [`trunc`] function.
    fn trunc(x: Self) -> Self;

    /// See the [`ceil`] function.
    fn ceil(x: Self) -> Self;

    /// See the [`floor`] function.
    fn floor(x: Self) -> Self;

    /// See the [`scalbn`] function.
    fn scalbn(x: Self, y: i32) -> Self;

    /// See the [`frexp`] function.
    fn frexp(x: Self) -> (Self, i32);

    /// See the [`hypot`] function.
    fn hypot(x: Self, y: Self) -> Self;

    /// See the [`sqrt`] function.
    fn sqrt(x: Self) -> Self;

    /// See the [`cbrt`] function.
    fn cbrt(x: Self) -> Self;

    /// See the [`exp`] function.
    fn exp(x: Self) -> Self;

    /// See the [`exp_m1`] function.
    fn exp_m1(x: Self) -> Self;

    /// See the [`exp2`] function.
    fn exp2(x: Self) -> Self;

    /// See the [`exp10`] function.
    fn exp10(x: Self) -> Self;

    /// See the [`log`] function.
    fn log(x: Self) -> Self;

    /// See the [`log_1p`] function.
    fn log_1p(x: Self) -> Self;

    /// See the [`log2`] function.
    fn log2(x: Self) -> Self;

    /// See the [`log10`] function.
    fn log10(x: Self) -> Self;

    /// See the [`pow`] function.
    fn pow(x: Self, y: Self) -> Self;

    /// See the [`powi`] function.
    fn powi(x: Self, y: i32) -> Self;

    /// See the [`sin`] function.
    fn sin(x: Self) -> Self;

    /// See the [`cos`] function.
    fn cos(x: Self) -> Self;

    /// See the [`sin_cos`] function.
    fn sin_cos(x: Self) -> (Self, Self);

    /// See the [`tan`] function.
    fn tan(x: Self) -> Self;

    /// See the [`sind`] function.
    fn sind(x: Self) -> Self;

    /// See the [`cosd`] function.
    fn cosd(x: Self) -> Self;

    /// See the [`sind_cosd`] function.
    fn sind_cosd(x: Self) -> (Self, Self);

    /// See the [`tand`] function.
    fn tand(x: Self) -> Self;

    /// See the [`sinpi`] function.
    fn sinpi(x: Self) -> Self;

    /// See the [`cospi`] function.
    fn cospi(x: Self) -> Self;

    /// See the [`sinpi_cospi`] function.
    fn sinpi_cospi(x: Self) -> (Self, Self);

    /// See the [`tanpi`] function.
    fn tanpi(x: Self) -> Self;

    /// See the [`asin`] function.
    fn asin(x: Self) -> Self;

    /// See the [`acos`] function.
    fn acos(x: Self) -> Self;

    /// See the [`atan`] function.
    fn atan(x: Self) -> Self;

    /// See the [`atan2`] function.
    fn atan2(y: Self, x: Self) -> Self;

    /// See the [`asind`] function.
    fn asind(x: Self) -> Self;

    /// See the [`acosd`] function.
    fn acosd(x: Self) -> Self;

    /// See the [`atand`] function.
    fn atand(x: Self) -> Self;

    /// See the [`atan2d`] function.
    fn atan2d(y: Self, x: Self) -> Self;

    /// See the [`asinpi`] function.
    fn asinpi(x: Self) -> Self;

    /// See the [`acospi`] function.
    fn acospi(x: Self) -> Self;

    /// See the [`atanpi`] function.
    fn atanpi(x: Self) -> Self;

    /// See the [`atan2pi`] function.
    fn atan2pi(y: Self, x: Self) -> Self;

    /// See the [`sinh`] function.
    fn sinh(x: Self) -> Self;

    /// See the [`cosh`] function.
    fn cosh(x: Self) -> Self;

    /// See the [`sinh_cosh`] function.
    fn sinh_cosh(x: Self) -> (Self, Self);

    /// See the [`tanh`] function.
    fn tanh(x: Self) -> Self;

    /// See the [`asinh`] function.
    fn asinh(x: Self) -> Self;

    /// See the [`acosh`] function.
    fn acosh(x: Self) -> Self;

    /// See the [`atanh`] function.
    fn atanh(x: Self) -> Self;

    /// See the [`tgamma`] function.
    fn tgamma(x: Self) -> Self;

    /// See the [`lgamma`] function.
    fn lgamma(x: Self) -> (Self, i8);
}

/// Calculates the absolute value of `x`
#[inline]
pub fn abs<F: FloatMath>(x: F) -> F {
    F::abs(x)
}

/// Returns a value with the magnitude of `x` and the sign of `y`
#[inline]
pub fn copysign<F: FloatMath>(x: F, y: F) -> F {
    F::copysign(x, y)
}

/// Rounds `x` to the nearest integer, ties round away from zero
pub fn round<F: FloatMath>(x: F) -> F {
    F::round(x)
}

/// Rounds `x` to the nearest integer that is not greater in magnitude than `x`
pub fn trunc<F: FloatMath>(x: F) -> F {
    F::trunc(x)
}

/// Rounds `x` to the nearest integer that is not less than `x`
pub fn ceil<F: FloatMath>(x: F) -> F {
    F::ceil(x)
}

/// Rounds `x` to the nearest integer that is not greater than `x`
pub fn floor<F: FloatMath>(x: F) -> F {
    F::floor(x)
}

/// Calculates `x` times two raised to `y`.
pub fn scalbn<F: FloatMath>(x: F, y: i32) -> F {
    F::scalbn(x, y)
}

/// Splits `x` into mantissa and exponent.
///
/// Returns `(m, e)` such as:
/// * `0.5 <= m < 1.0`
/// * `x = m * 2^e`
///
/// When `x` is zero, infinity or NaN, returns `x` as mantissa and zero as
/// exponent.
pub fn frexp<F: FloatMath>(x: F) -> (F, i32) {
    F::frexp(x)
}

/// Calculates the Pythagorean addition of `x` and `y` with and error of less
/// than 1 ULP
///
/// The Pythagorean addition of `x` and `y` is equal to the length of the
/// hypotenuse of a triangle with sides of length `x` and `y`.
///
/// Special cases:
/// * Returns positive infinity if `x` or `y` is infinity
/// * Returns NaN if `x` or `y` is NaN and neither is infinity
pub fn hypot<F: FloatMath>(x: F, y: F) -> F {
    F::hypot(x, y)
}

/// Calculates the square root of `x` with an error of less than 0.5 ULP.
///
/// Special cases:
/// * Returns negative zero if `x` is negative zero
/// * Returns positive infinity if `x` is positive infinity
/// * Returns NaN if `x` is NaN or negative non-zero (including infinity)
pub fn sqrt<F: FloatMath>(x: F) -> F {
    F::sqrt(x)
}

/// Calculates the cube root of `x` with and error of less than 1 ULP.
///
/// Special cases:
/// * Returns negative zero if `x` is negative zero
/// * Returns positive infinity if `x` is positive infinity
/// * Returns negative infinity if `x` is negative infinity
/// * Returns NaN if `x` is NaN
pub fn cbrt<F: FloatMath>(x: F) -> F {
    F::cbrt(x)
}

/// Calculates Euler's number raised to `x` with an error of less than 1 ULP
///
/// Special cases:
/// * Returns positive infinity if `x` is positive infinity
/// * Returns zero if `x` is negative infinity
/// * Returns NaN if `x` is NaN
pub fn exp<F: FloatMath>(x: F) -> F {
    F::exp(x)
}

/// Calculates `exp(x) - 1.0` with an error of less than 1 ULP
///
/// Special cases:
/// * Returns negative zero if `x` is negative zero
/// * Returns positive infinity if `x` is positive infinity
/// * Returns minus one if `x` is negative infinity
/// * Returns NaN if `x` is NaN
pub fn exp_m1<F: FloatMath>(x: F) -> F {
    F::exp_m1(x)
}

/// Calculates 2 raised to `x` with an error of less than 1 ULP
///
/// Special cases:
/// * Returns positive infinity if `x` is positive infinity
/// * Returns zero if `x` is negative infinity
/// * Returns NaN if `x` is NaN
pub fn exp2<F: FloatMath>(x: F) -> F {
    F::exp2(x)
}

/// Calculates 10 raised to `x` with an error of less than 1 ULP
///
/// Special cases:
/// * Returns positive infinity if `x` is positive infinity
/// * Returns zero if `x` is negative infinity
/// * Returns NaN if `x` is NaN
pub fn exp10<F: FloatMath>(x: F) -> F {
    F::exp10(x)
}

/// Calculates the natural logarithm of `x` with an error of less than 1 ULP
///
/// Special cases:
/// * Returns negative infinity if `x` is positive or negative zero
/// * Returns positive infinity if `x` is positive infinity
/// * Returns NaN if `x` is NaN or negative non-zero (including infinity)
pub fn log<F: FloatMath>(x: F) -> F {
    F::log(x)
}

/// Calculates the natural logarithm of `x + 1` with an error of less than 1 ULP
///
/// Special cases:
/// * Returns negative zero if `x` is negative zero
/// * Returns negative infinity if `x` is minus one
/// * Returns positive infinity if `x` is positive infinity
/// * Returns NaN if `x` is NaN or less than minus one (including negative
///   infinity)
pub fn log_1p<F: FloatMath>(x: F) -> F {
    F::log_1p(x)
}

/// Calculates the base-2 logarithm of `x` with an error of less than 1 ULP
///
/// Special cases:
/// * Returns negative infinity if `x` is positive or negative zero
/// * Returns positive infinity if `x` is positive infinity
/// * Returns NaN if `x` is NaN or negative non-zero (including infinity)
pub fn log2<F: FloatMath>(x: F) -> F {
    F::log2(x)
}

/// Calculates the base-10 logarithm of `x` with an error of less than 1 ULP
///
/// Special cases:
/// * Returns negative infinity if `x` is positive or negative zero
/// * Returns positive infinity if `x` is positive infinity
/// * Returns NaN if `x` is NaN or negative non-zero (including infinity)
pub fn log10<F: FloatMath>(x: F) -> F {
    F::log10(x)
}

/// Calculates `x` raised to `y` with an error of less than 1 ULP
///
/// Special cases:
/// * Returns 1 when `x` is 1 or `y` is zero
/// * Returns 1 when `x` is -1 and `y` is positive or negative infinity
/// * Returns NaN when `x` is NaN and `y` is not zero
/// * Returns NaN when `y` is NaN and `z` is not 1
/// * Returns positive zero when `x` is positive zero and `y` is positive
/// * Returns positive zero when `x` is positive infinity and `y` is negative
/// * Returns positive infinity when `x` is positive infinity and `y` is
///   positive
/// * Returns positive infinity when `x` is positive zero and `y` is negative
/// * Returns positive zero when `x` is negative zero and `y` is positive and
///   not an odd integer
/// * Returns positive zero when `x` is negative infinity and `y` is negative
///   and not an odd integer
/// * Returns positive infinity when `x` is negative infinity and `y` is
///   positive and not an odd integer
/// * Returns positive infinity when `x` is negative zero and `y` is negative
///   and not an odd integer
/// * Returns negative zero when `x` is negative zero and `y` is positive and an
///   odd integer
/// * Returns negative zero when `x` is negative infinity and `y` is negative
///   and an odd integer
/// * Returns negative infinity when `x` is negative infinity and `y` is
///   positive and an odd integer
/// * Returns negative infinity when `x` is negative zero and `y` is negative
///   and an odd integer
/// * Returns positive zero when the absolute value of `x` is less than one and
///   `y` is positive infinity
/// * Returns positive zero when the absolute value of `x` is greater than one
///   and `y` is negative infinity
/// * Returns positive infinity when the absolute value of `x` is less than one
///   and `y` is negative infinity
/// * Returns positive infinity when the absolute value of `x` is greater than
///   one and `y` is positive infinity
/// * Returns NaN when `x` is negative and `y` is finite and not integer
pub fn pow<F: FloatMath>(x: F, y: F) -> F {
    F::pow(x, y)
}

/// Calculates `x` raised to `y` with an error of less than 1 ULP
///
/// Special cases:
/// * Returns 1 when `x` is 1 or `y` is zero
/// * Returns NaN when `x` is NaN and `y` is not zero
/// * Returns NaN when `y` is NaN and `z` is not 1
/// * Returns positive zero when `x` is positive zero and `y` is positive
/// * Returns positive zero when `x` is positive infinity and `y` is negative
/// * Returns positive infinity when `x` is positive infinity and `y` is
///   positive
/// * Returns positive infinity when `x` is positive zero and `y` is negative
/// * Returns positive zero when `x` is negative zero and `y` is positive and
///   even
/// * Returns positive zero when `x` is negative infinity and `y` is negative
///   and even
/// * Returns positive infinity when `x` is negative infinity and `y` is
///   positive and even
/// * Returns positive infinity when `x` is negative zero and `y` is negative
///   and even
/// * Returns negative zero when `x` is negative zero and `y` is positive and
///   odd
/// * Returns negative zero when `x` is negative infinity and `y` is negative
///   and odd
/// * Returns negative infinity when `x` is negative infinity and `y` is
///   positive and odd
/// * Returns negative infinity when `x` is negative zero and `y` is negative
///   and odd
pub fn powi<F: FloatMath>(x: F, y: i32) -> F {
    F::powi(x, y)
}

/// Calculates the sine of `x` radians with an error of less than 1 ULP
///
/// Special cases:
/// * Returns negative zero if `x` is negative zero
/// * Returns NaN if `x` is infinity or NaN
pub fn sin<F: FloatMath>(x: F) -> F {
    F::sin(x)
}

/// Calculates the cosine of `x` radians with an error of less than 1 ULP
///
/// Special cases:
/// * Returns NaN if `x` is infinity or NaN
pub fn cos<F: FloatMath>(x: F) -> F {
    F::cos(x)
}

/// Calculates the sine and the cosine of `x` radians
///
/// The same accuracy and special cases of [`sin`] and [`cos`] also
/// apply to this function. Using this function can be faster than
/// using [`sin`] and [`cos`] separately.
pub fn sin_cos<F: FloatMath>(x: F) -> (F, F) {
    F::sin_cos(x)
}

/// Calculates the tangent of `x` radians with an error of less than 1 ULP
///
/// Special cases:
/// * Returns negative zero if `x` is negative zero
/// * Returns NaN if `x` is infinity or NaN
pub fn tan<F: FloatMath>(x: F) -> F {
    F::tan(x)
}

/// Calculates the sine of `x` degrees with an error of less than 1 ULP
///
/// Special cases:
/// * Returns negative zero if `x` is negative zero
/// * Returns NaN if `x` is infinity or NaN
pub fn sind<F: FloatMath>(x: F) -> F {
    F::sind(x)
}

/// Calculates the cosine of `x` degrees with an error of less than 1 ULP
///
/// Special cases:
/// * Returns NaN if `x` is infinity or NaN
pub fn cosd<F: FloatMath>(x: F) -> F {
    F::cosd(x)
}

/// Calculates the sine and the cosine of `x` degrees
///
/// The same accuracy and special cases of [`sind`] and [`cosd`] also apply to
/// this function. Using this function can be faster than using [`sind`] and
/// [`cosd`] separately.
pub fn sind_cosd<F: FloatMath>(x: F) -> (F, F) {
    F::sind_cosd(x)
}

/// Calculates the tangent of `x` degrees with an error of less than 1 ULP
///
/// Special cases:
/// * Returns negative zero if `x` is negative zero
/// * Returns NaN if `x` is infinity or NaN
pub fn tand<F: FloatMath>(x: F) -> F {
    F::tand(x)
}

/// Calculates the sine of `x` half-revolutions with an error of less
/// than 1 ULP
///
/// Special cases:
/// * Returns negative zero if `x` is negative zero
/// * Returns NaN if `x` is infinity or NaN
pub fn sinpi<F: FloatMath>(x: F) -> F {
    F::sinpi(x)
}

/// Calculates the cosine of `x` half-revolutions with an error of
/// less than 1 ULP
///
/// Special cases:
/// * Returns NaN if `x` is infinity or NaN
pub fn cospi<F: FloatMath>(x: F) -> F {
    F::cospi(x)
}

/// Calculates the sine and the cosine of `x` half-revolutions
///
/// The same accuracy and special cases of [`sinpi`] and [`cospi`] also apply to
/// this function. Using this function can be faster than using [`sinpi`] and
/// [`cospi`] separately.
pub fn sinpi_cospi<F: FloatMath>(x: F) -> (F, F) {
    F::sinpi_cospi(x)
}

/// Calculates the tangent of `x` half-revolutions with an error of less than 1
/// ULP
///
/// Special cases:
/// * Returns negative zero if `x` is negative zero
/// * Returns NaN if `x` is infinity or NaN
pub fn tanpi<F: FloatMath>(x: F) -> F {
    F::tanpi(x)
}

/// Calculates the arcsine of `x`, returning the result in radians, with an
/// error of less than 1 ULP
///
/// Special cases:
/// * Returns negative zero if `x` is negative zero
/// * Returns NaN if `x` is NaN or greater than one in magnitude (including
///   infinity)
pub fn asin<F: FloatMath>(x: F) -> F {
    F::asin(x)
}

/// Calculates the arccosine of `x`, returning the result in radians, with an
/// error of less than 1 ULP
///
/// Special cases:
/// * Returns NaN if `x` is NaN or greater than one in magnitude (including
///   infinity)
pub fn acos<F: FloatMath>(x: F) -> F {
    F::acos(x)
}

/// Calculates the arctangent of `x`, returning the result in radians,
/// with an error of less than 1 ULP
///
/// Special cases:
/// * Returns negative zero if `x` is negative zero
/// * Returns NaN if `x` is NaN
/// * Returns π/2 if `x` is positive infinity
/// * Returns -π/2 if `x` is negative infinity
pub fn atan<F: FloatMath>(x: F) -> F {
    F::atan(x)
}

/// Calculates the 2-argument arctangent of `x` and 'y', returning the
/// result in radians with an error of less than 1 ULP
///
/// Special cases:
/// * Returns NaN if `x` is NaN or `y` is NaN
/// * Returns positive zero if `y` is positive zero `x` is positive (zero,
///   finite or infinity)
/// * Returns negative zero if `y` is negative zero `x` is positive (zero,
///   finite or infinity)
/// * Returns π if `y` is positive zero `x` is positive (zero, finite or
///   infinity)
/// * Returns -π if `y` is negative zero `x` is negative (zero, finite or
///   infinity)
/// * Returns π/2 if `y` is positive infinity and `x` is zero or finite
/// * Returns -π/2 if `y` is negative infinity and `x` is zero or finite
/// * Returns positive zero if `x` is positive infinity and `y` is positive
///   (zero or finite)
/// * Returns negative zero if `x` is positive infinity and `y` is negative
///   (zero or finite)
/// * Returns π if `x` is negative infinity and `y` is positive (zero or finite)
/// * Returns -π if `x` is negative infinity and `y` is negative (zero or
///   finite)
/// * Returns π/4 if `x` is positive infinity and `y` is positive infinity
/// * Returns -π/4 if `x` is positive infinity and `y` is negative infinity
/// * Returns 3π/4 if `x` is negative infinity and `y` is positive infinity
/// * Returns -3π/4 if `x` is negative infinity and `y` is negative infinity
pub fn atan2<F: FloatMath>(y: F, x: F) -> F {
    F::atan2(y, x)
}

/// Calculates the arcsine of `x`, returning the result in degrees, with an
/// error of less than 1 ULP
///
/// Special cases:
/// * Returns negative zero if `x` is negative zero
/// * Returns NaN if `x` is NaN or greater than one in magnitude (including
///   infinity)
pub fn asind<F: FloatMath>(x: F) -> F {
    F::asind(x)
}

/// Calculates the arccosine of `x`, returning the result in degrees, with an
/// error of less than 1 ULP
///
/// Special cases:
/// * Returns NaN if `x` is NaN or greater than one in magnitude (including
///   infinity)
pub fn acosd<F: FloatMath>(x: F) -> F {
    F::acosd(x)
}

/// Calculates the arctangent of `x`, returning the result in degrees, with an
/// error of less than 1 ULP
///
/// Special cases:
/// * Returns negative zero if `x` is negative zero
/// * Returns NaN if `x` is NaN
/// * Returns 90 if `x` is positive infinity
/// * Returns -90 if `x` is negative infinity
pub fn atand<F: FloatMath>(x: F) -> F {
    F::atand(x)
}

/// Calculates the 2-argument arctangent of `x` and 'y', returning the result in
/// degrees, with an error of less than 1 ULP
///
/// Special cases:
/// * Returns NaN if `x` is NaN or `y` is NaN
/// * Returns positive zero if `y` is positive zero `x` is positive (zero,
///   finite or infinity)
/// * Returns negative zero if `y` is negative zero `x` is positive (zero,
///   finite or infinity)
/// * Returns 180 if `y` is positive zero `x` is positive (zero, finite or
///   infinity)
/// * Returns -180 if `y` is negative zero `x` is negative (zero, finite or
///   infinity)
/// * Returns 90 if `y` is positive infinity and `x` is zero or finite
/// * Returns -90 if `y` is negative infinity and `x` is zero or finite
/// * Returns positive zero if `x` is positive infinity and `y` is positive
///   (zero or finite)
/// * Returns negative zero if `x` is positive infinity and `y` is negative
///   (zero or finite)
/// * Returns 180 if `x` is negative infinity and `y` is positive (zero or
///   finite)
/// * Returns -180 if `x` is negative infinity and `y` is negative (zero or
///   finite)
/// * Returns 45 if `x` is positive infinity and `y` is positive infinity
/// * Returns -45 if `x` is positive infinity and `y` is negative infinity
/// * Returns 135 if `x` is negative infinity and `y` is positive infinity
/// * Returns -135 if `x` is negative infinity and `y` is negative infinity
pub fn atan2d<F: FloatMath>(y: F, x: F) -> F {
    F::atan2d(y, x)
}

/// Calculates the arcsine of `x`, returning the result in half-revolutions,
/// with an error of less than 1 ULP
///
/// Special cases:
/// * Returns negative zero if `x` is negative zero
/// * Returns NaN if `x` is NaN or greater than one in magnitude (including
///   infinity)
pub fn asinpi<F: FloatMath>(x: F) -> F {
    F::asinpi(x)
}

/// Calculates the arccosine of `x`, returning the result in half-revolutions,
/// with an error of less than 1 ULP
///
/// Special cases:
/// * Returns NaN if `x` is NaN or greater than one in magnitude (including
///   infinity)
pub fn acospi<F: FloatMath>(x: F) -> F {
    F::acospi(x)
}

/// Calculates the arctangent of `x`, returning the result in half-revolutions,
/// with an error of less than 1 ULP
///
/// Special cases:
/// * Returns negative zero if `x` is negative zero
/// * Returns NaN if `x` is NaN
/// * Returns 0.5 if `x` is positive infinity
/// * Returns -0.5 if `x` is negative infinity
pub fn atanpi<F: FloatMath>(x: F) -> F {
    F::atanpi(x)
}

/// Calculates the 2-argument arctangent of `x` and 'y', returning the result in
/// half-revolutions, with an error of less than 1 ULP
///
/// Special cases:
/// * Returns NaN if `x` is NaN or `y` is NaN
/// * Returns positive zero if `y` is positive zero `x` is positive (zero,
///   finite or infinity)
/// * Returns negative zero if `y` is negative zero `x` is positive (zero,
///   finite or infinity)
/// * Returns 1 if `y` is positive zero `x` is positive (zero, finite or
///   infinity)
/// * Returns -1 if `y` is negative zero `x` is negative (zero, finite or
///   infinity)
/// * Returns 0.5 if `y` is positive infinity and `x` is zero or finite
/// * Returns -0.5 if `y` is negative infinity and `x` is zero or finite
/// * Returns positive zero if `x` is positive infinity and `y` is positive
///   (zero or finite)
/// * Returns negative zero if `x` is positive infinity and `y` is negative
///   (zero or finite)
/// * Returns 1 if `x` is negative infinity and `y` is positive (zero or finite)
/// * Returns -1 if `x` is negative infinity and `y` is negative (zero or
///   finite)
/// * Returns 0.25 if `x` is positive infinity and `y` is positive infinity
/// * Returns -0.25 if `x` is positive infinity and `y` is negative infinity
/// * Returns 0.75 if `x` is negative infinity and `y` is positive infinity
/// * Returns -0.75 if `x` is negative infinity and `y` is negative infinity
pub fn atan2pi<F: FloatMath>(y: F, x: F) -> F {
    F::atan2pi(y, x)
}

/// Calculates the hyperbolic sine of `x` with an error of less than 1 ULP
///
/// Special cases:
/// * Returns negative zero if `x` is negative zero
/// * Returns positive infinity if `x` is positive infinity
/// * Returns negative infinity if `x` is negative infinity
/// * Returns NaN if `x` is NaN
pub fn sinh<F: FloatMath>(x: F) -> F {
    F::sinh(x)
}

/// Calculates the hyperbolic cosine of `x` with an error of less than 1 ULP
///
/// Special cases:
/// * Returns positive infinity if `x` is positive or negative infinity
/// * Returns NaN if `x` is NaN
pub fn cosh<F: FloatMath>(x: F) -> F {
    F::cosh(x)
}

/// Calculates the hyerbolic sine and hyerbolic cosine of `x`
///
/// The same accuracy and special cases of [`sinh`] and [`cosh`] also apply to
/// this function. Using this function can be faster than using [`sinh`] and
/// [`cosh`] separately.
pub fn sinh_cosh<F: FloatMath>(x: F) -> (F, F) {
    F::sinh_cosh(x)
}

/// Calculates the hyperbolic tangent of `x` with an error of less than 1 ULP
///
/// Special cases:
/// * Returns negative zero if `x` is negative zero
/// * Returns one if `x` is positive infinity
/// * Returns minus one if `x` is negative infinity
/// * Returns NaN if `x` is NaN
pub fn tanh<F: FloatMath>(x: F) -> F {
    F::tanh(x)
}

/// Calculates the hyperbolic arcsine of `x` with an error of less than 1 ULP
///
/// Special cases:
/// * Returns NaN if `x` is NaN
/// * Returns negative zero if `x` is negative zero
/// * Returns positive infinity if `x` is positive infinity
/// * Returns negative infinity if `x` is negative infinity
pub fn asinh<F: FloatMath>(x: F) -> F {
    F::asinh(x)
}

/// Calculates the hyperbolic arccosine of `x` with an error of less than 1 ULP
///
/// Special cases:
/// * Returns NaN if `x` is NaN or less than one (including negative infinity)
/// * Returns positive infinity if `x` is positive infinity
pub fn acosh<F: FloatMath>(x: F) -> F {
    F::acosh(x)
}

/// Calculates the hyperbolic arctangent of `x` with an error of less than 1 ULP
///
/// Special cases:
/// * Returns NaN if `x` is NaN or greater than 1 in magnitude
/// * Returns negative zero if `x` is negative zero
/// * Returns positive infinity if `x` is 1
/// * Returns negative infinity if `x` is -1
pub fn atanh<F: FloatMath>(x: F) -> F {
    F::atanh(x)
}

/// Calculates the gamma function of `x`
///
/// When `x` is greater than 0.5, the error is less than 1 ULP, otherwise the
/// error is less than 2 ULP.
///
/// Special cases:
/// * Returns NaN if `x` is NaN, negative infinity or a negative integer
/// * Returns positive infinity if `x` is positive infinity
/// * Returns positive infinity if `x` is positive zero
/// * Returns negative infinity if `x` is negative zero
pub fn tgamma<F: FloatMath>(x: F) -> F {
    F::tgamma(x)
}

/// Calculates the logarithm of the absolute value of the gamma function of `x`
///
/// The integer field of the returned tuple is `1` when the gamma function of
/// `x` is positive, `-1` when the gamma function of `x` is negative, and `0`
/// when the sign of the gamma function of `x` is not defined.
///
/// The error is less than 2 ULP in most cases. However, for some negative
/// values of `x`, the error can be in the order of 200 ULP.
///
/// Special cases:
/// * Returns NaN if `x` is NaN or negative infinity
/// * Returns positive infinity if `x` is positive infinity
/// * Returns positive infinity if `x` is zero or a negative integer
///
/// The sign is considered undefined when `x` is NaN, a negative integer or
/// negative infinity.
pub fn lgamma<F: FloatMath>(x: F) -> (F, i8) {
    F::lgamma(x)
}
