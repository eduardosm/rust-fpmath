//! Inverse hyperbolic functions for `f32`.
//!
//! They use the same formulas as the `f64` functions (see
//! `crate::f64::inv_hyperbolic`), but evaluated with plain `f64` arithmetic,
//! and with the `ln` of `f32` (see `super::log`).
//!
//! `x^2` does not overflow in `f64`, but `asinh` and `acosh` still use
//! `ln(2 * x) + ln((1 + sqrt(1 ± 1 / x^2)) / 2)` when `x >= 16`, which avoids
//! the square root.

use super::log::ln_f64;
use crate::f64::{acosh_large_corr, asinh_large_corr, fast_sqrt};
use crate::traits::Float as _;

impl crate::generic::InvHyperbolic for f32 {
    #[inline]
    fn asinh_finite(x: Self) -> Self {
        let x_exp = x.exponent();
        if x_exp < -12 {
            // asinh(x) = x - x^3 / 6 + ..., where |x| < 2^-12, so
            // x^3 / 6 - ... < 2^-26 * |x|, which is less than half the spacing
            // of `f32` around x, and the result rounds to x
            x
        } else if x_exp < -1 {
            // |x| < 0.5
            asinh_small(f64::from(x)) as f32
        } else {
            // asinh(x) = asinh(|x|) * sgn(x)
            let ax = f64::from(x.abs());
            let r = if x_exp >= 4 {
                // |x| >= 16
                ln_f64(ax, 1) + asinh_large_corr(1.0 / (ax * ax))
            } else {
                // `x^2 + 1` and the sum are rounded, the square root has the
                // error of `fast_sqrt`, and `ln_f64` adds its absolute error
                // for arguments of up to 53 bits, while the result is greater
                // than 0.48
                ln_f64(ax + fast_sqrt(ax * ax + 1.0).0, 0)
            };
            (r as f32).copysign(x)
        }
    }

    #[inline]
    fn acosh_finite(x: Self) -> Self {
        let x = f64::from(x);
        let r = if x >= 16.0 {
            ln_f64(x, 1) + acosh_large_corr(1.0 / (x * x))
        } else {
            // `x^2 - 1` is exact (`x^2` has at most 48 bits), the square
            // root has the error of `fast_sqrt`, the sum is rounded, and
            // `ln_f64` adds its absolute error for arguments of up to 53
            // bits, while the result is greater than 2^-12
            // (`x >= 1 + 2^-23`)
            ln_f64(x + fast_sqrt(x * x - 1.0).0, 0)
        };
        r as f32
    }

    #[inline]
    fn atanh_finite(x: Self) -> Self {
        let x_exp = x.exponent();
        if x_exp < -12 {
            // atanh(x) = x + x^3 / 3 + ..., where |x| < 2^-12, so
            // x^3 / 3 + ... < 2^-25.5 * |x|, which is less than half the
            // spacing of `f32` around x, and the result rounds to x
            x
        } else if x_exp < -5 {
            // |x| < 2^-5
            atanh_small(f64::from(x)) as f32
        } else {
            // atanh(x) = ln((1 + |x|) / (1 - |x|)) / 2 * sgn(x)
            // `1 + |x|` and `1 - |x|` are exact, the division is rounded, and
            // `ln_f64` adds its absolute error for arguments of up to 53
            // bits, while the result is greater than 2^-5
            let ax = f64::from(x.abs());
            let r = 0.5 * ln_f64((1.0 + ax) / (1.0 - ax), 0);
            (r as f32).copysign(x)
        }
    }
}

/// Returns `asinh(x)`, for `|x| < 0.5`, with a relative error of about
/// 2^-38.5.
#[inline]
fn asinh_small(x: f64) -> f64 {
    // asinh(x) ~= x + x^3 * (K3 + K5 * x^2 + ... + K15 * x^12)
    // GENERATE: asinh_poly f64 7 1e-30 0.5
    const K3: f64 = f64::from_bits(0xBFC5555553955811); // -1.6666666585177883e-1
    const K5: f64 = f64::from_bits(0x3FB33331C49DC1A6); // 7.49999146480663e-2
    const K7: f64 = f64::from_bits(0xBFA6DB07D7FB972B); // -4.463982116230477e-2
    const K9: f64 = f64::from_bits(0x3F9F0EE998C23DC2); // 3.033032412400672e-2
    const K11: f64 = f64::from_bits(0xBF966BB6A96E3296); // -2.1895269480909883e-2
    const K13: f64 = f64::from_bits(0x3F8E6A0EA4D62209); // 1.4850725560578775e-2
    const K15: f64 = f64::from_bits(0xBF7AFF0872EE6213); // -6.590874675474828e-3

    let x2 = x * x;
    let x4 = x2 * x2;
    let x8 = x4 * x4;
    let q = ((K3 + x2 * K5) + x4 * (K7 + x2 * K9)) + x8 * ((K11 + x2 * K13) + x4 * K15);
    x + x * (x2 * q)
}

/// Returns `atanh(x)`, for `|x| < 2^-5`, with a relative error of about
/// 2^-50.
#[inline]
fn atanh_small(x: f64) -> f64 {
    // atanh(x) ~= x + x^3 * (K3 + K5 * x^2 + K7 * x^4)
    // GENERATE: atanh_poly f64 3 -0.0313 0.0313
    const K3: f64 = f64::from_bits(0x3FD55555555A9CAB); // 3.3333333335253695e-1
    const K5: f64 = f64::from_bits(0x3FC999989A355109); // 1.999998810739572e-1
    const K7: f64 = f64::from_bits(0x3FC250008F1163FF); // 1.4306647287119742e-1

    let x2 = x * x;
    x + x * (x2 * (K3 + x2 * (K5 + x2 * K7)))
}
