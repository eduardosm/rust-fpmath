//! Error functions for `f32`.
//!
//! They use the same formulas as the `f64` functions (see `crate::f64::erf`),
//! but evaluated with plain `f64` arithmetic:
//!
//! * `erf(x) = x * P(x^2)` for `|x| < POLY_MAX`, with a shorter polynomial.
//! * `erfc(|x|) = exp(-x^2) * erfcx(|x|)` otherwise, where `x^2` is exact,
//!   the exponential is the one of `super::exp` and `erfcx` uses the
//!   polynomials of the `f64` functions without their last terms (see
//!   `crate::f64::erfcx_short`).
//!
//! `erfcx` also uses those polynomials up to `ASYMPT_MIN`, and then a shorter
//! polynomial for the asymptotic series (see `erfcx_asympt`).

use super::exp::exp_f64;
use crate::f64::erfcx_short;
use crate::traits::Float as _;

/// `erf(x)` is calculated as `x * P(x^2)` for `|x| < POLY_MAX`, and from
/// `erfcx(|x|)` otherwise.
const POLY_MAX: f64 = 0.5;

/// `erfc(x) < 2^-25.9` for `x >= SATURATION_MIN`, so `erf(x)` rounds to one
/// and `erfc(-x)` rounds to two.
const SATURATION_MIN: f64 = 4.0;

/// `erfc(x) < 2^-155` for `x >= ERFC_ZERO_MIN`, so it rounds to zero.
const ERFC_ZERO_MIN: f64 = 10.25;

/// `erfcx(x)` is calculated with the polynomials of `erfcx_short` for
/// `POLY_MAX <= x < ASYMPT_MIN`, and as `Q(1 / x^2) / x` otherwise.
const ASYMPT_MIN: f32 = 28.0;

/// `2 * exp(x^2) > 2^128` for `x > ~9.39`, so `erfcx(-x)` overflows for
/// `x >= ERFCX_INF_MIN`.
const ERFCX_INF_MIN: f32 = 9.5;

impl crate::generic::Erf for f32 {
    #[inline]
    fn erf_finite(x: Self) -> Self {
        // erf(x) = erf(|x|) * sgn(x)
        let ax = f64::from(x.abs());
        let r = if ax < POLY_MAX {
            erf_small(ax)
        } else if ax < SATURATION_MIN {
            1.0 - erfc_large(ax)
        } else {
            1.0
        };
        (r as f32).copysign(x)
    }

    #[inline]
    fn erfc_finite(x: Self) -> Self {
        let x = f64::from(x);
        let ax = x.abs();
        let r = if ax < POLY_MAX {
            // erfc(x) = 1 - erf(x), where |erf(x)| < 0.53
            1.0 - erf_small(x)
        } else if x.is_sign_negative() {
            if ax < SATURATION_MIN {
                // erfc(x) = 2 - erfc(|x|)
                2.0 - erfc_large(ax)
            } else {
                2.0
            }
        } else if x < ERFC_ZERO_MIN {
            // It can be subnormal as `f32`, but it is normal as `f64`
            erfc_large(x)
        } else {
            0.0
        };
        r as f32
    }

    #[inline]
    fn erfcx_finite(x: Self) -> Self {
        // The comparisons use `f32` values and each branch converts its
        // argument to `f64`: converting `x` first is much slower in a loop on
        // x86-64, where the conversion would depend on the previous contents
        // of its destination register.
        let ax = x.abs();
        let r = if ax < POLY_MAX as f32 {
            // erfcx(x) = exp(x^2) * (1 - erf(x)), where |erf(x)| < 0.53 and
            // `x^2` is exact
            let x = f64::from(x);
            exp_f64(x * x) * (1.0 - erf_small(x))
        } else if x.is_sign_negative() {
            if ax < ERFCX_INF_MIN {
                // erfcx(x) = 2 * exp(x^2) - erfcx(|x|), where `x^2` is exact
                // and erfcx(|x|) < 0.62. It can overflow as `f32`, but not as
                // `f64`.
                let ax = f64::from(ax);
                2.0 * exp_f64(ax * ax) - erfcx_short(ax)
            } else {
                f64::INFINITY
            }
        } else if x < ASYMPT_MIN {
            erfcx_short(f64::from(x))
        } else {
            // It can be subnormal as `f32`, but it is normal as `f64`
            erfcx_asympt(f64::from(x))
        };
        r as f32
    }
}

/// Returns `erf(x)`, for `|x| < POLY_MAX`, with a relative error of about
/// 2^-43.
#[inline]
fn erf_small(x: f64) -> f64 {
    // erf(x) ~= x * (K0 + K1 * x^2 + ... + K6 * x^12)
    // GENERATE: erf_poly f64 7 0.5
    const K0: f64 = f64::from_bits(0x3FF20DD75042999A); // 1.1283791670954089e0
    const K1: f64 = f64::from_bits(0xBFD812746AF85FCB); // -3.7612638899144873e-1
    const K2: f64 = f64::from_bits(0x3FBCE2F20EF51CAA); // 1.1283791413471769e-1
    const K3: f64 = f64::from_bits(0xBF9B82CA0CD320C5); // -2.686610892558929e-2
    const K4: f64 = f64::from_bits(0x3F7564FF20BB8D75); // 5.223270981989149e-3
    const K5: f64 = f64::from_bits(0xBF4BDFDA05D01DB5); // -8.506598056791577e-4
    const K6: f64 = f64::from_bits(0x3F1C5D8DB19EC58B); // 1.0820557909487474e-4

    if x.exponent() < -27 {
        // |x| < 2^-27, so the terms after `K0` are less than 2^-55 * K0 and
        // do not change it. Skipping them also avoids the powers of `x` that
        // are subnormal (in `f64`) when `x` is a subnormal `f32`, which can
        // be very slow on some CPUs.
        return x * K0;
    }

    let x2 = x * x;
    let x4 = x2 * x2;
    x * (((K0 + x2 * K1) + x4 * (K2 + x2 * K3)) + (x4 * x4) * ((K4 + x2 * K5) + x4 * K6))
}

/// Returns `erfc(x)`, for `POLY_MAX <= x < 11.3` (so `x^2 < 128`), where `x`
/// has at most 26 significant bits (so `x^2` is exact), with a relative error
/// of about 2^-39.5:
/// * `exp(-x^2)`: 2^-40.5 (see `exp_f64`).
/// * `erfcx(x)`: 2^-40.5 (see `erfcx_short`).
#[inline]
fn erfc_large(x: f64) -> f64 {
    exp_f64(-(x * x)) * erfcx_short(x)
}

/// Returns `erfcx(x)`, for `x >= ASYMPT_MIN`, with a relative error of about
/// 2^-42.5.
#[inline]
fn erfcx_asympt(x: f64) -> f64 {
    // erfcx(x) ~= (K0 + K1 / x^2 + K2 / x^4 + K3 / x^6) / x
    // GENERATE: erfcx_asympt_poly f64 4 28
    const K0: f64 = f64::from_bits(0x3FE20DD750429B6D); // 5.641895835477563e-1
    const K1: f64 = f64::from_bits(0xBFD20DD74EC3DA32); // -2.8209479038142427e-1
    const K2: f64 = f64::from_bits(0x3FDB14A723A98BC8); // 4.231355522408049e-1
    const K3: f64 = f64::from_bits(0xBFF0C81D3D19E20D); // -1.0488560091745314e0

    let t = 1.0 / x;
    let s = t * t;
    t * ((K0 + s * K1) + (s * s) * (K2 + s * K3))
}
