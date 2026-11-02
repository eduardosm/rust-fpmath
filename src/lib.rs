#![warn(
    rust_2018_idioms,
    trivial_casts,
    trivial_numeric_casts,
    unreachable_pub,
    unused_qualifications
)]
#![no_std]

//! A pure-Rust floating point math library.
//!
//! This crate is `no_std`.
//!
//! The following math functions are implemented:
//!
//! * Sign operations ([`abs`], [`copysign`]).
//! * Rounding ([`round`], [`trunc`], [`ceil`], [`floor`]).
//! * Exponential ([`exp`], [`exp_m1`], [`exp2`], [`exp2_m1`], [`exp10`],
//!   [`exp10_m1`]).
//! * Logarithmic ([`ln`], [`ln_1p`], [`log2`], [`log2_1p`], [`log10`],
//!   [`log10_1p`]).
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
//! * Gamma ([`gamma`], [`ln_gamma`]).
//! * Error function ([`erf`], [`erfc`], [`erfcx`]).
//! * Bessel functions
//!   - First kind ([`j0`], [`j1`], [`jn`]).
//!   - Second kind ([`y0`], [`y1`], [`yn`]).
//!
//! All functions are implemted for the native floating point types [`prim@f32`]
//! and [`prim@f64`].
//!
//! The [`FloatMath`] trait is used to identify types that support the math
//! functions.
//!
//! # Accuracy
//!
//! The results of most functions have an error of less than 0.5 + 2<sup>-10</sup>
//! ULP (units in the last place), with these exceptions:
//!
//! * [`sqrt`], whose result is always correctly rounded (its error is less
//!   than 0.5 ULP).
//! * [`jn`] and [`yn`] with an order `n` other than 0, 1 and -1, when the
//!   absolute value of `x` is greater than the one of `n`. There, the
//!   functions oscillate around zero with an amplitude of
//!   `sqrt(jn(n, x)^2 + yn(n, x)^2)`, and their error can be larger by up to
//!   2<sup>-63</sup> times that amplitude, which is only noticeable in the
//!   results that are much smaller than it (close to a zero of the function).

// Uncomment to use `dbg!`
//extern crate std;

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

mod f32;
mod f64;
mod generic;
mod int;
mod traits;

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

    /// See the [`exp2_m1`] function.
    fn exp2_m1(x: Self) -> Self;

    /// See the [`exp10`] function.
    fn exp10(x: Self) -> Self;

    /// See the [`exp10_m1`] function.
    fn exp10_m1(x: Self) -> Self;

    /// See the [`ln`] function.
    fn ln(x: Self) -> Self;

    /// See the [`ln_1p`] function.
    fn ln_1p(x: Self) -> Self;

    /// See the [`log2`] function.
    fn log2(x: Self) -> Self;

    /// See the [`log2_1p`] function.
    fn log2_1p(x: Self) -> Self;

    /// See the [`log10`] function.
    fn log10(x: Self) -> Self;

    /// See the [`log10_1p`] function.
    fn log10_1p(x: Self) -> Self;

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

    /// See the [`gamma`] function.
    fn gamma(x: Self) -> Self;

    /// See the [`ln_gamma`] function.
    fn ln_gamma(x: Self) -> (Self, i8);

    /// See the [`erf`] function.
    fn erf(x: Self) -> Self;

    /// See the [`erfc`] function.
    fn erfc(x: Self) -> Self;

    /// See the [`erfcx`] function.
    fn erfcx(x: Self) -> Self;

    /// See the [`j0`] function.
    fn j0(x: Self) -> Self;

    /// See the [`j1`] function.
    fn j1(x: Self) -> Self;

    /// See the [`jn`] function.
    fn jn(n: i32, x: Self) -> Self;

    /// See the [`y0`] function.
    fn y0(x: Self) -> Self;

    /// See the [`y1`] function.
    fn y1(x: Self) -> Self;

    /// See the [`yn`] function.
    fn yn(n: i32, x: Self) -> Self;
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
/// * `0.5 <= |m| < 1.0`
/// * `x = m * 2^e`
///
/// When `x` is zero, infinity or NaN, returns `x` as mantissa and zero as
/// exponent.
pub fn frexp<F: FloatMath>(x: F) -> (F, i32) {
    F::frexp(x)
}

/// Calculates the Pythagorean addition of `x` and `y`
///
/// The Pythagorean addition of `x` and `y` is equal to the length of the
/// hypotenuse of a triangle with sides of length `x` and `y`.
///
/// Special cases:
///
/// | `x`                 | `y`                 | Result |
/// |---------------------|---------------------|--------|
/// | ±∞                  | any (including NaN) | +∞     |
/// | any (including NaN) | ±∞                  | +∞     |
/// | NaN                 | not ±∞              | NaN    |
/// | not ±∞              | NaN                 | NaN    |
pub fn hypot<F: FloatMath>(x: F, y: F) -> F {
    F::hypot(x, y)
}

/// Calculates the square root of `x` with an error of less than 0.5 ULP.
///
/// Special cases:
///
/// | `x`                | Result |
/// |--------------------|--------|
/// | -0                 | -0     |
/// | +∞                 | +∞     |
/// | < 0 (including -∞) | NaN    |
/// | NaN                | NaN    |
pub fn sqrt<F: FloatMath>(x: F) -> F {
    F::sqrt(x)
}

/// Calculates the cube root of `x`.
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | -0  | -0     |
/// | +∞  | +∞     |
/// | -∞  | -∞     |
/// | NaN | NaN    |
pub fn cbrt<F: FloatMath>(x: F) -> F {
    F::cbrt(x)
}

/// Calculates Euler's number raised to `x`
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | +∞  | +∞     |
/// | -∞  | +0     |
/// | NaN | NaN    |
pub fn exp<F: FloatMath>(x: F) -> F {
    F::exp(x)
}

/// Calculates `exp(x) - 1.0`
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | -0  | -0     |
/// | +∞  | +∞     |
/// | -∞  | -1     |
/// | NaN | NaN    |
pub fn exp_m1<F: FloatMath>(x: F) -> F {
    F::exp_m1(x)
}

/// Calculates 2 raised to `x`
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | +∞  | +∞     |
/// | -∞  | +0     |
/// | NaN | NaN    |
pub fn exp2<F: FloatMath>(x: F) -> F {
    F::exp2(x)
}

/// Calculates `exp2(x) - 1.0`
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | -0  | -0     |
/// | +∞  | +∞     |
/// | -∞  | -1     |
/// | NaN | NaN    |
pub fn exp2_m1<F: FloatMath>(x: F) -> F {
    F::exp2_m1(x)
}

/// Calculates 10 raised to `x`
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | +∞  | +∞     |
/// | -∞  | +0     |
/// | NaN | NaN    |
pub fn exp10<F: FloatMath>(x: F) -> F {
    F::exp10(x)
}

/// Calculates `exp10(x) - 1.0`
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | -0  | -0     |
/// | +∞  | +∞     |
/// | -∞  | -1     |
/// | NaN | NaN    |
pub fn exp10_m1<F: FloatMath>(x: F) -> F {
    F::exp10_m1(x)
}

/// Calculates the natural logarithm of `x`
///
/// Special cases:
///
/// | `x`                | Result |
/// |--------------------|--------|
/// | ±0                 | -∞     |
/// | +∞                 | +∞     |
/// | < 0 (including -∞) | NaN    |
/// | NaN                | NaN    |
pub fn ln<F: FloatMath>(x: F) -> F {
    F::ln(x)
}

/// Calculates the natural logarithm of `x + 1`
///
/// Special cases:
///
/// | `x`                 | Result |
/// |---------------------|--------|
/// | -0                  | -0     |
/// | -1                  | -∞     |
/// | +∞                  | +∞     |
/// | < -1 (including -∞) | NaN    |
/// | NaN                 | NaN    |
pub fn ln_1p<F: FloatMath>(x: F) -> F {
    F::ln_1p(x)
}

/// Calculates the base-2 logarithm of `x`
///
/// Special cases:
///
/// | `x`                | Result |
/// |--------------------|--------|
/// | ±0                 | -∞     |
/// | +∞                 | +∞     |
/// | < 0 (including -∞) | NaN    |
/// | NaN                | NaN    |
pub fn log2<F: FloatMath>(x: F) -> F {
    F::log2(x)
}

/// Calculates the base-2 logarithm of `x + 1`
///
/// Special cases:
///
/// | `x`                 | Result |
/// |---------------------|--------|
/// | -0                  | -0     |
/// | -1                  | -∞     |
/// | +∞                  | +∞     |
/// | < -1 (including -∞) | NaN    |
/// | NaN                 | NaN    |
pub fn log2_1p<F: FloatMath>(x: F) -> F {
    F::log2_1p(x)
}

/// Calculates the base-10 logarithm of `x`
///
/// Special cases:
///
/// | `x`                | Result |
/// |--------------------|--------|
/// | ±0                 | -∞     |
/// | +∞                 | +∞     |
/// | < 0 (including -∞) | NaN    |
/// | NaN                | NaN    |
pub fn log10<F: FloatMath>(x: F) -> F {
    F::log10(x)
}

/// Calculates the base-10 logarithm of `x + 1`
///
/// Special cases:
///
/// | `x`                 | Result |
/// |---------------------|--------|
/// | -0                  | -0     |
/// | -1                  | -∞     |
/// | +∞                  | +∞     |
/// | < -1 (including -∞) | NaN    |
/// | NaN                 | NaN    |
pub fn log10_1p<F: FloatMath>(x: F) -> F {
    F::log10_1p(x)
}

/// Calculates `x` raised to `y`
///
/// Special cases:
///
/// | `x`                 | `y`                     | Result |
/// |---------------------|-------------------------|--------|
/// | 1                   | any (including NaN)     | 1      |
/// | any (including NaN) | ±0                      | 1      |
/// | -1                  | ±∞                      | 1      |
/// | NaN                 | not ±0                  | NaN    |
/// | not 1               | NaN                     | NaN    |
/// | +0                  | > 0                     | +0     |
/// | +0                  | < 0                     | +∞     |
/// | -0                  | > 0, odd integer        | -0     |
/// | -0                  | > 0, not an odd integer | +0     |
/// | -0                  | < 0, odd integer        | -∞     |
/// | -0                  | < 0, not an odd integer | +∞     |
/// | +∞                  | > 0                     | +∞     |
/// | +∞                  | < 0                     | +0     |
/// | -∞                  | > 0, odd integer        | -∞     |
/// | -∞                  | > 0, not an odd integer | +∞     |
/// | -∞                  | < 0, odd integer        | -0     |
/// | -∞                  | < 0, not an odd integer | +0     |
/// | > -1 and < 1        | +∞                      | +0     |
/// | > -1 and < 1        | -∞                      | +∞     |
/// | < -1 or > 1         | +∞                      | +∞     |
/// | < -1 or > 1         | -∞                      | +0     |
/// | < 0, finite         | finite, not an integer  | NaN    |
pub fn pow<F: FloatMath>(x: F, y: F) -> F {
    F::pow(x, y)
}

/// Calculates `x` raised to `y`
///
/// Special cases:
///
/// | `x`                 | `y`       | Result |
/// |---------------------|-----------|--------|
/// | 1                   | any       | 1      |
/// | any (including NaN) | 0         | 1      |
/// | NaN                 | not 0     | NaN    |
/// | +0                  | > 0       | +0     |
/// | +0                  | < 0       | +∞     |
/// | -0                  | > 0, odd  | -0     |
/// | -0                  | > 0, even | +0     |
/// | -0                  | < 0, odd  | -∞     |
/// | -0                  | < 0, even | +∞     |
/// | +∞                  | > 0       | +∞     |
/// | +∞                  | < 0       | +0     |
/// | -∞                  | > 0, odd  | -∞     |
/// | -∞                  | > 0, even | +∞     |
/// | -∞                  | < 0, odd  | -0     |
/// | -∞                  | < 0, even | +0     |
pub fn powi<F: FloatMath>(x: F, y: i32) -> F {
    F::powi(x, y)
}

/// Calculates the sine of `x` radians
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | -0  | -0     |
/// | ±∞  | NaN    |
/// | NaN | NaN    |
pub fn sin<F: FloatMath>(x: F) -> F {
    F::sin(x)
}

/// Calculates the cosine of `x` radians
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | ±∞  | NaN    |
/// | NaN | NaN    |
pub fn cos<F: FloatMath>(x: F) -> F {
    F::cos(x)
}

/// Calculates the sine and the cosine of `x` radians
///
/// The same special cases of [`sin`] and [`cos`] also apply to this function.
/// Using this function can be faster than using [`sin`] and [`cos`] separately.
pub fn sin_cos<F: FloatMath>(x: F) -> (F, F) {
    F::sin_cos(x)
}

/// Calculates the tangent of `x` radians
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | -0  | -0     |
/// | ±∞  | NaN    |
/// | NaN | NaN    |
pub fn tan<F: FloatMath>(x: F) -> F {
    F::tan(x)
}

/// Calculates the sine of `x` degrees
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | -0  | -0     |
/// | ±∞  | NaN    |
/// | NaN | NaN    |
pub fn sind<F: FloatMath>(x: F) -> F {
    F::sind(x)
}

/// Calculates the cosine of `x` degrees
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | ±∞  | NaN    |
/// | NaN | NaN    |
pub fn cosd<F: FloatMath>(x: F) -> F {
    F::cosd(x)
}

/// Calculates the sine and the cosine of `x` degrees
///
/// The same special cases of [`sind`] and [`cosd`] also apply to this function.
/// Using this function can be faster than using [`sind`] and [`cosd`]
/// separately.
pub fn sind_cosd<F: FloatMath>(x: F) -> (F, F) {
    F::sind_cosd(x)
}

/// Calculates the tangent of `x` degrees
///
/// Special cases:
///
/// | `x`                    | Result |
/// |------------------------|--------|
/// | -0                     | -0     |
/// | 90 + 360k (integer k)  | +∞     |
/// | -90 + 360k (integer k) | -∞     |
/// | ±∞                     | NaN    |
/// | NaN                    | NaN    |
pub fn tand<F: FloatMath>(x: F) -> F {
    F::tand(x)
}

/// Calculates the sine of `x` half-revolutions
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | -0  | -0     |
/// | ±∞  | NaN    |
/// | NaN | NaN    |
pub fn sinpi<F: FloatMath>(x: F) -> F {
    F::sinpi(x)
}

/// Calculates the cosine of `x` half-revolutions
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | ±∞  | NaN    |
/// | NaN | NaN    |
pub fn cospi<F: FloatMath>(x: F) -> F {
    F::cospi(x)
}

/// Calculates the sine and the cosine of `x` half-revolutions
///
/// The same special cases of [`sinpi`] and [`cospi`] also apply to this
/// function. Using this function can be faster than using [`sinpi`] and
/// [`cospi`] separately.
pub fn sinpi_cospi<F: FloatMath>(x: F) -> (F, F) {
    F::sinpi_cospi(x)
}

/// Calculates the tangent of `x` half-revolutions
///
/// Special cases:
///
/// | `x`                   | Result |
/// |-----------------------|--------|
/// | -0                    | -0     |
/// | 0.5 + 2k (integer k)  | +∞     |
/// | -0.5 + 2k (integer k) | -∞     |
/// | ±∞                    | NaN    |
/// | NaN                   | NaN    |
pub fn tanpi<F: FloatMath>(x: F) -> F {
    F::tanpi(x)
}

/// Calculates the arcsine of `x`, returning the result in radians
///
/// Special cases:
///
/// | `x`                        | Result |
/// |----------------------------|--------|
/// | -0                         | -0     |
/// | < -1 or > 1 (including ±∞) | NaN    |
/// | NaN                        | NaN    |
pub fn asin<F: FloatMath>(x: F) -> F {
    F::asin(x)
}

/// Calculates the arccosine of `x`, returning the result in radians
///
/// Special cases:
///
/// | `x`                        | Result |
/// |----------------------------|--------|
/// | < -1 or > 1 (including ±∞) | NaN    |
/// | NaN                        | NaN    |
pub fn acos<F: FloatMath>(x: F) -> F {
    F::acos(x)
}

/// Calculates the arctangent of `x`, returning the result in radians
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | -0  | -0     |
/// | +∞  | π/2    |
/// | -∞  | -π/2   |
/// | NaN | NaN    |
pub fn atan<F: FloatMath>(x: F) -> F {
    F::atan(x)
}

/// Calculates the 2-argument arctangent of `y` and `x`, returning the result in
/// radians
///
/// Special cases:
///
/// | `y`         | `x`                      | Result |
/// |-------------|--------------------------|--------|
/// | NaN         | any                      | NaN    |
/// | any         | NaN                      | NaN    |
/// | +0          | +0 or > 0 (including +∞) | +0     |
/// | -0          | +0 or > 0 (including +∞) | -0     |
/// | +0          | -0 or < 0 (including -∞) | π      |
/// | -0          | -0 or < 0 (including -∞) | -π     |
/// | > 0, finite | +∞                       | +0     |
/// | < 0, finite | +∞                       | -0     |
/// | > 0, finite | -∞                       | π      |
/// | < 0, finite | -∞                       | -π     |
/// | > 0, finite | ±0                       | π/2    |
/// | < 0, finite | ±0                       | -π/2   |
/// | +∞          | finite (including ±0)    | π/2    |
/// | -∞          | finite (including ±0)    | -π/2   |
/// | +∞          | +∞                       | π/4    |
/// | -∞          | +∞                       | -π/4   |
/// | +∞          | -∞                       | 3π/4   |
/// | -∞          | -∞                       | -3π/4  |
pub fn atan2<F: FloatMath>(y: F, x: F) -> F {
    F::atan2(y, x)
}

/// Calculates the arcsine of `x`, returning the result in degrees
///
/// Special cases:
///
/// | `x`                        | Result |
/// |----------------------------|--------|
/// | -0                         | -0     |
/// | < -1 or > 1 (including ±∞) | NaN    |
/// | NaN                        | NaN    |
pub fn asind<F: FloatMath>(x: F) -> F {
    F::asind(x)
}

/// Calculates the arccosine of `x`, returning the result in degrees
///
/// Special cases:
///
/// | `x`                        | Result |
/// |----------------------------|--------|
/// | < -1 or > 1 (including ±∞) | NaN    |
/// | NaN                        | NaN    |
pub fn acosd<F: FloatMath>(x: F) -> F {
    F::acosd(x)
}

/// Calculates the arctangent of `x`, returning the result in degrees
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | -0  | -0     |
/// | +∞  | 90     |
/// | -∞  | -90    |
/// | NaN | NaN    |
pub fn atand<F: FloatMath>(x: F) -> F {
    F::atand(x)
}

/// Calculates the 2-argument arctangent of `y` and `x`, returning the result in
/// degrees
///
/// Special cases:
///
/// | `y`         | `x`                      | Result |
/// |-------------|--------------------------|--------|
/// | NaN         | any                      | NaN    |
/// | any         | NaN                      | NaN    |
/// | +0          | +0 or > 0 (including +∞) | +0     |
/// | -0          | +0 or > 0 (including +∞) | -0     |
/// | +0          | -0 or < 0 (including -∞) | 180    |
/// | -0          | -0 or < 0 (including -∞) | -180   |
/// | > 0, finite | +∞                       | +0     |
/// | < 0, finite | +∞                       | -0     |
/// | > 0, finite | -∞                       | 180    |
/// | < 0, finite | -∞                       | -180   |
/// | > 0, finite | ±0                       | 90     |
/// | < 0, finite | ±0                       | -90    |
/// | +∞          | finite (including ±0)    | 90     |
/// | -∞          | finite (including ±0)    | -90    |
/// | +∞          | +∞                       | 45     |
/// | -∞          | +∞                       | -45    |
/// | +∞          | -∞                       | 135    |
/// | -∞          | -∞                       | -135   |
pub fn atan2d<F: FloatMath>(y: F, x: F) -> F {
    F::atan2d(y, x)
}

/// Calculates the arcsine of `x`, returning the result in half-revolutions
///
/// Special cases:
///
/// | `x`                        | Result |
/// |----------------------------|--------|
/// | -0                         | -0     |
/// | < -1 or > 1 (including ±∞) | NaN    |
/// | NaN                        | NaN    |
pub fn asinpi<F: FloatMath>(x: F) -> F {
    F::asinpi(x)
}

/// Calculates the arccosine of `x`, returning the result in half-revolutions
///
/// Special cases:
///
/// | `x`                        | Result |
/// |----------------------------|--------|
/// | < -1 or > 1 (including ±∞) | NaN    |
/// | NaN                        | NaN    |
pub fn acospi<F: FloatMath>(x: F) -> F {
    F::acospi(x)
}

/// Calculates the arctangent of `x`, returning the result in half-revolutions
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | -0  | -0     |
/// | +∞  | 0.5    |
/// | -∞  | -0.5   |
/// | NaN | NaN    |
pub fn atanpi<F: FloatMath>(x: F) -> F {
    F::atanpi(x)
}

/// Calculates the 2-argument arctangent of `y` and `x`, returning the result in
/// half-revolutions
///
/// Special cases:
///
/// | `y`         | `x`                      | Result |
/// |-------------|--------------------------|--------|
/// | NaN         | any                      | NaN    |
/// | any         | NaN                      | NaN    |
/// | +0          | +0 or > 0 (including +∞) | +0     |
/// | -0          | +0 or > 0 (including +∞) | -0     |
/// | +0          | -0 or < 0 (including -∞) | 1      |
/// | -0          | -0 or < 0 (including -∞) | -1     |
/// | > 0, finite | +∞                       | +0     |
/// | < 0, finite | +∞                       | -0     |
/// | > 0, finite | -∞                       | 1      |
/// | < 0, finite | -∞                       | -1     |
/// | > 0, finite | ±0                       | 0.5    |
/// | < 0, finite | ±0                       | -0.5   |
/// | +∞          | finite (including ±0)    | 0.5    |
/// | -∞          | finite (including ±0)    | -0.5   |
/// | +∞          | +∞                       | 0.25   |
/// | -∞          | +∞                       | -0.25  |
/// | +∞          | -∞                       | 0.75   |
/// | -∞          | -∞                       | -0.75  |
pub fn atan2pi<F: FloatMath>(y: F, x: F) -> F {
    F::atan2pi(y, x)
}

/// Calculates the hyperbolic sine of `x`
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | -0  | -0     |
/// | +∞  | +∞     |
/// | -∞  | -∞     |
/// | NaN | NaN    |
pub fn sinh<F: FloatMath>(x: F) -> F {
    F::sinh(x)
}

/// Calculates the hyperbolic cosine of `x`
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | ±∞  | +∞     |
/// | NaN | NaN    |
pub fn cosh<F: FloatMath>(x: F) -> F {
    F::cosh(x)
}

/// Calculates the hyperbolic sine and hyperbolic cosine of `x`
///
/// The same special cases of [`sinh`] and [`cosh`] also apply to this function.
/// Using this function can be faster than using [`sinh`] and [`cosh`]
/// separately.
pub fn sinh_cosh<F: FloatMath>(x: F) -> (F, F) {
    F::sinh_cosh(x)
}

/// Calculates the hyperbolic tangent of `x`
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | -0  | -0     |
/// | +∞  | 1      |
/// | -∞  | -1     |
/// | NaN | NaN    |
pub fn tanh<F: FloatMath>(x: F) -> F {
    F::tanh(x)
}

/// Calculates the hyperbolic arcsine of `x`
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | -0  | -0     |
/// | +∞  | +∞     |
/// | -∞  | -∞     |
/// | NaN | NaN    |
pub fn asinh<F: FloatMath>(x: F) -> F {
    F::asinh(x)
}

/// Calculates the hyperbolic arccosine of `x`
///
/// Special cases:
///
/// | `x`                | Result |
/// |--------------------|--------|
/// | +∞                 | +∞     |
/// | < 1 (including -∞) | NaN    |
/// | NaN                | NaN    |
pub fn acosh<F: FloatMath>(x: F) -> F {
    F::acosh(x)
}

/// Calculates the hyperbolic arctangent of `x`
///
/// Special cases:
///
/// | `x`                        | Result |
/// |----------------------------|--------|
/// | -0                         | -0     |
/// | 1                          | +∞     |
/// | -1                         | -∞     |
/// | < -1 or > 1 (including ±∞) | NaN    |
/// | NaN                        | NaN    |
pub fn atanh<F: FloatMath>(x: F) -> F {
    F::atanh(x)
}

/// Calculates the gamma function of `x`
///
/// When `x` is negative and not an integer, the sign of the gamma function
/// alternates between consecutive integers, and results that underflow to zero
/// keep that sign. For example, `gamma(-180.5)` returns negative zero and
/// `gamma(-181.5)` returns positive zero.
///
/// Special cases:
///
/// | `x`          | Result |
/// |--------------|--------|
/// | +0           | +∞     |
/// | -0           | -∞     |
/// | +∞           | +∞     |
/// | -∞           | NaN    |
/// | < 0, integer | NaN    |
/// | NaN          | NaN    |
pub fn gamma<F: FloatMath>(x: F) -> F {
    F::gamma(x)
}

/// Calculates the logarithm of the absolute value of the gamma function of `x`
///
/// The integer field of the returned tuple is `1` when the gamma function of
/// `x` is positive, `-1` when the gamma function of `x` is negative, and `0`
/// when the sign of the gamma function of `x` is not defined.
///
/// Special cases:
///
/// | `x`          | Result | Sign |
/// |--------------|--------|------|
/// | +0           | +∞     | 1    |
/// | -0           | +∞     | -1   |
/// | +∞           | +∞     | 1    |
/// | -∞           | NaN    | 0    |
/// | < 0, integer | +∞     | 0    |
/// | NaN          | NaN    | 0    |
pub fn ln_gamma<F: FloatMath>(x: F) -> (F, i8) {
    F::ln_gamma(x)
}

/// Calculates the error function of `x`
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | -0  | -0     |
/// | +∞  | 1      |
/// | -∞  | -1     |
/// | NaN | NaN    |
pub fn erf<F: FloatMath>(x: F) -> F {
    F::erf(x)
}

/// Calculates the complementary error function of `x`
///
/// The complementary error function of `x` is `1 - erf(x)`, which this
/// function calculates without losing accuracy when `erf(x)` is close to one.
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | +∞  | +0     |
/// | -∞  | 2      |
/// | NaN | NaN    |
pub fn erfc<F: FloatMath>(x: F) -> F {
    F::erfc(x)
}

/// Calculates the scaled complementary error function of `x`
///
/// The scaled complementary error function of `x` is `exp(x^2) * erfc(x)`,
/// which this function calculates without the underflow of `erfc(x)` and the
/// overflow of `exp(x^2)` when `x` is large.
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | +∞  | +0     |
/// | -∞  | +∞     |
/// | NaN | NaN    |
pub fn erfcx<F: FloatMath>(x: F) -> F {
    F::erfcx(x)
}

/// Calculates the Bessel function of the first kind of order 0 of `x`
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | ±∞  | +0     |
/// | NaN | NaN    |
pub fn j0<F: FloatMath>(x: F) -> F {
    F::j0(x)
}

/// Calculates the Bessel function of the first kind of order 1 of `x`
///
/// Special cases:
///
/// | `x` | Result |
/// |-----|--------|
/// | -0  | -0     |
/// | +∞  | +0     |
/// | -∞  | -0     |
/// | NaN | NaN    |
pub fn j1<F: FloatMath>(x: F) -> F {
    F::j1(x)
}

/// Calculates the Bessel function of the first kind of order `n` of `x`
///
/// `jn(-n, x)` and `jn(n, -x)` are equal to `jn(n, x)` when `n` is even and to
/// `-jn(n, x)` when `n` is odd, including the sign of the zero results.
///
/// When the absolute value of `n` is greater than one, the results close to
/// zero can be less accurate than the ones of other functions (see the
/// [accuracy](crate#accuracy) of the crate).
///
/// Special cases:
///
/// | `n`         | `x`      | Result |
/// |-------------|----------|--------|
/// | 0           | ±0       | 1      |
/// | even, not 0 | ±0       | +0     |
/// | even        | ±∞       | +0     |
/// | odd, > 0    | +0 or +∞ | +0     |
/// | odd, > 0    | -0 or -∞ | -0     |
/// | odd, < 0    | +0 or +∞ | -0     |
/// | odd, < 0    | -0 or -∞ | +0     |
/// | any         | NaN      | NaN    |
pub fn jn<F: FloatMath>(n: i32, x: F) -> F {
    F::jn(n, x)
}

/// Calculates the Bessel function of the second kind of order 0 of `x`
///
/// Special cases:
///
/// | `x`                | Result |
/// |--------------------|--------|
/// | ±0                 | -∞     |
/// | +∞                 | +0     |
/// | < 0 (including -∞) | NaN    |
/// | NaN                | NaN    |
pub fn y0<F: FloatMath>(x: F) -> F {
    F::y0(x)
}

/// Calculates the Bessel function of the second kind of order 1 of `x`
///
/// Special cases:
///
/// | `x`                | Result |
/// |--------------------|--------|
/// | ±0                 | -∞     |
/// | +∞                 | +0     |
/// | < 0 (including -∞) | NaN    |
/// | NaN                | NaN    |
pub fn y1<F: FloatMath>(x: F) -> F {
    F::y1(x)
}

/// Calculates the Bessel function of the second kind of order `n` of `x`
///
/// `yn(-n, x)` is equal to `yn(n, x)` when `n` is even and to `-yn(n, x)` when
/// `n` is odd.
///
/// When the absolute value of `n` is greater than one, the results close to
/// zero can be less accurate than the ones of other functions (see the
/// [accuracy](crate#accuracy) of the crate).
///
/// Special cases:
///
/// | `n`         | `x`                | Result |
/// |-------------|--------------------|--------|
/// | even or > 0 | ±0                 | -∞     |
/// | odd, < 0    | ±0                 | +∞     |
/// | even or > 0 | +∞                 | +0     |
/// | odd, < 0    | +∞                 | -0     |
/// | any         | < 0 (including -∞) | NaN    |
/// | any         | NaN                | NaN    |
pub fn yn<F: FloatMath>(n: i32, x: F) -> F {
    F::yn(n, x)
}
