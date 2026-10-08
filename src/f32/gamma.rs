//! Gamma functions for `f32`.
//!
//! They use the same formulas as the `f64` functions (see
//! `crate::f64::gamma`), but evaluated with plain `f64` arithmetic, with
//! shorter polynomials. The factors of the products of consecutive values are
//! exact (except `x + 1` and `x + 2`, which can be rounded when `|x| < 2^-28`),
//! so each multiplication (or rounded factor) adds a relative error of at
//! most 2^-53.
//!
//! `gamma` also uses the Stirling series for large arguments (with the
//! exponential of `super::exp`) and the reflection formula for very negative
//! arguments, but evaluated as `π / (sin(π * x) * Γ(1 - x))`.
//!
//! `ln_gamma` calculates its results close to zero (around its zeros with
//! negative arguments) like the `f64` function (see
//! `crate::f64::ln_gamma_near_zero`).

use super::exp::exp_f64;
use super::log::ln_f64_accurate;
use super::trigonometric::sinpi_f64;
use crate::f64::{ln_gamma_near_zero, round_i32};
use crate::traits::Float as _;

// GENERATE: consts f64 PI LN_PI
const PI: f64 = f64::from_bits(0x400921FB54442D18); // 3.141592653589793e0
const LN_PI: f64 = f64::from_bits(0x3FF250D048E7A1BD); // 1.1447298858494002e0

/// `gamma` uses the Stirling series for `x >= GAMMA_STIRLING_MIN`.
const GAMMA_STIRLING_MIN: f64 = 10.0;

/// `gamma` uses the reflection formula for `x <= GAMMA_REFLECTION`.
const GAMMA_REFLECTION: f64 = -10.0;

/// `ln_gamma` uses the Stirling series for `x >= LN_GAMMA_STIRLING_MIN`.
const LN_GAMMA_STIRLING_MIN: f64 = 4.0;

/// `ln_gamma` uses the reflection formula for `x <= LN_GAMMA_REFLECTION`,
/// where `|ln|Γ(x)|| > 1.2`.
const LN_GAMMA_REFLECTION: f64 = -10.0;

impl crate::generic::Gamma for f32 {
    #[inline]
    fn gamma_finite(x: Self) -> Self {
        if x >= 36.0 {
            // Overflow
            return f32::INFINITY;
        } else if x <= -43.0 {
            // Underflow, with the sign of `Γ(x)`, which is negative between
            // `-2 * k - 1` and `-2 * k` (`x` is not an integer)
            let xi = (-x) as u32;
            return 0.0.set_sign((xi & 1) == 0);
        }

        let xd = f64::from(x);
        let r = if xd >= GAMMA_STIRLING_MIN {
            // ln(Γ(x)) < 92.2 has an absolute error of about 2^-44.5, and
            // the exponential adds a relative error of about 2^-40.5
            exp_f64(stirling(xd, -0.5))
        } else if xd > GAMMA_REFLECTION {
            gamma_small(xd)
        } else {
            // Γ(x) = π / (sin(π * x) * Γ(1 - x)), where
            // Γ(1 - x) = exp(ln(Γ(1 + y))), with y = -x, is less than Γ(44)
            PI / (sinpi_f64(x) * exp_f64(stirling(-xd, 0.5)))
        };
        r as f32
    }

    #[inline]
    fn ln_gamma_finite(x: Self) -> (Self, i8) {
        let xd = f64::from(x);
        if xd >= LN_GAMMA_STIRLING_MIN {
            // The absolute error is about 2^-46 (with a relative error of
            // about 2^-52 for large arguments), and the result is greater
            // than 1.79
            (stirling(xd, -0.5) as f32, 1)
        } else if xd > LN_GAMMA_REFLECTION {
            ln_gamma_small(xd)
        } else {
            // ln|Γ(x)| = ln(π) - ln|sin(π * x)| - ln(Γ(1 + y)), with y = -x,
            // where the result is greater than 1.2 in magnitude (the zeros of
            // `ln|Γ(x)|` are between the representable arguments), and the
            // absolute error is about 2^-45
            let s = sinpi_f64(x);
            let sign = if s.is_sign_negative() { -1 } else { 1 };
            let r = (LN_PI - ln_f64_accurate(s.abs(), 0)) - stirling(-xd, 0.5);
            (r as f32, sign)
        }
    }
}

/// Returns `Γ(x)`, for `GAMMA_REFLECTION < x < GAMMA_STIRLING_MIN`, where `x`
/// is not zero or a negative integer, with a relative error of about 2^-42.5.
#[inline]
fn gamma_small(x: f64) -> f64 {
    // x = n + f, where `f` is exact
    let (nf, n) = round_i32(x);
    let f = x - nf;
    let g = gamma_poly(f);

    if n >= 3 {
        // Γ(x) = Γ(3 + f) * (3 + f) * ... * (x - 1)
        let mut p = g;
        let mut v = x;
        for _ in 3..n {
            v -= 1.0;
            p *= v;
        }
        p
    } else {
        // Γ(x) = Γ(3 + f) / (x * (x + 1) * ... * (2 + f))
        let mut p = x;
        let mut v = x;
        for _ in n..2 {
            v += 1.0;
            p *= v;
        }
        g / p
    }
}

/// Returns `(ln|Γ(x)|, sign(Γ(x)))`, for
/// `LN_GAMMA_REFLECTION < x < LN_GAMMA_STIRLING_MIN`, where `x` is not zero or
/// a negative integer.
#[inline]
fn ln_gamma_small(x: f64) -> (f32, i8) {
    // x = n + f, where `f` is exact
    let (nf, n) = round_i32(x);
    let f = x - nf;
    let g = ln_gamma_poly(f);

    if n == 2 {
        (g as f32, 1)
    } else if n > 2 {
        // ln(Γ(x)) = ln(Γ(2 + f)) + ln((2 + f) * ... * (x - 1)), where the
        // product is greater than 1.5 and has at most 48 bits (so the
        // logarithm has a relative error of about 2^-52), and the result is
        // greater than 0.28
        let mut v = x - 1.0;
        let mut p = v;
        for _ in 3..n {
            v -= 1.0;
            p *= v;
        }
        ((g + ln_f64_accurate(p, 0)) as f32, 1)
    } else {
        // ln|Γ(x)| = ln(Γ(2 + f)) - ln|x * (x + 1) * ... * (1 + f)|
        let mut p = x;
        let mut v = x;
        for _ in n..1 {
            v += 1.0;
            p *= v;
        }
        let sign = if p.is_sign_negative() { -1 } else { 1 };
        let r = g - ln_f64_accurate(p.abs(), 0);
        // When `n <= -2`, the absolute error is about 2^-45.8 (mostly from
        // `ln(Γ(2 + f))`), which is not small enough relative to results close
        // to zero (around the zeros between -10 and -2), so they are
        // calculated like in the `f64` function for such results. When
        // `n = 1`, the result (close to zero when `x` is close to 1) is about
        // `-γ * f`, and both terms have small relative errors.
        let r = if n <= -2 && r.abs() < f64::exp2i_fast(-12) {
            ln_gamma_near_zero(x)
        } else {
            r
        };
        (r as f32, sign)
    }
}

/// Returns `ln(Γ(y + c + 1/2))`, where `c = ±1/2` (so it is `ln(Γ(y))` or
/// `ln(Γ(1 + y))`), for `y >= 4`, with the Stirling series:
///
/// `ln(Γ(y + c + 1/2)) = (y + c) * ln(y) - y + K0 + K1 / y + K3 / y^3 + ...`
///
/// The absolute error is about 2^-46 (mostly from the polynomial), plus about
/// 2^-52 relative to `(y + c) * ln(y)` (from the logarithm and the rounding of
/// the product).
#[inline]
fn stirling(y: f64, c: f64) -> f64 {
    // Coefficients fitted for `y >= 4`, with an absolute error of 2^-46.1
    // GENERATE: gamma_stirling_poly f64 7 4
    const K0: f64 = f64::from_bits(0x3FED67F1C864BEB5); // 9.189385332046728e-1
    const K1: f64 = f64::from_bits(0x3FB5555555555555); // 8.333333333333333e-2
    const K3: f64 = f64::from_bits(0xBF66C16C064F015F); // -2.777777658110533e-3
    const K5: f64 = f64::from_bits(0x3F4A016E868F5E2F); // 7.936277080444249e-4
    const K7: f64 = f64::from_bits(0xBF43740A11D30A4D); // -5.936669509981146e-4
    const K9: f64 = f64::from_bits(0x3F49F1F10EFAD7C3); // 7.917811435251765e-4
    const K11: f64 = f64::from_bits(0xBF5251C4E8BF406D); // -1.118128103360784e-3

    let l = ln_f64_accurate(y, 0);
    let t = 1.0 / y;
    let s = t * t;
    let s2 = s * s;
    let q = t * (K1 + s * ((K3 + s * K5) + s2 * ((K7 + s * K9) + s2 * K11)));
    // `y + c` is exact when `y < 2^52`, otherwise its rounding error is less
    // than 2^-53 relative to it
    ((y + c) * l - y) + (K0 + q)
}

/// Returns `Γ(3 + f)`, for `|f| <= 1/2`, with a relative error of about
/// 2^-42.5.
#[inline]
fn gamma_poly(f: f64) -> f64 {
    // GENERATE: gamma_poly f64 11 3 -0.5 0.5
    const K0: f64 = f64::from_bits(0x3FFFFFFFFFFFFB63); // 1.9999999999997378e0
    const K1: f64 = f64::from_bits(0x3FFD877303907B0B); // 1.8455686701997276e0
    const K2: f64 = f64::from_bits(0x3FF3F1854798EA11); // 1.2464649960298895e0
    const K3: f64 = f64::from_bits(0x3FE2665A2BBE04D6); // 5.749941686955491e-1
    const K4: f64 = f64::from_bits(0x3FCD731875001A06); // 2.3007493698987974e-1
    const K5: f64 = f64::from_bits(0x3FB2DEFD580FADB4); // 7.371505165035935e-2
    const K6: f64 = f64::from_bits(0x3F9691F63404ED60); // 2.2041174817088005e-2
    const K7: f64 = f64::from_bits(0x3F765162965081FE); // 5.448708635441777e-3
    const K8: f64 = f64::from_bits(0x3F56321078C16967); // 1.3547097069222748e-3
    const K9: f64 = f64::from_bits(0x3F315D4D52B8118B); // 2.64960649218681e-4
    const K10: f64 = f64::from_bits(0x3F1085DFAD4EC625); // 6.303003040277666e-5
    const K11: f64 = f64::from_bits(0x3EE16073DE7D1CA5); // 8.285888682414632e-6

    let f2 = f * f;
    let f4 = f2 * f2;
    let f8 = f4 * f4;
    ((K0 + f * K1) + f2 * (K2 + f * K3))
        + f4 * ((K4 + f * K5) + f2 * (K6 + f * K7))
        + f8 * ((K8 + f * K9) + f2 * (K10 + f * K11))
}

/// Returns `ln(Γ(2 + f))`, for `|f| <= 1/2`, with a relative error of about
/// 2^-44.
#[inline]
fn ln_gamma_poly(f: f64) -> f64 {
    // GENERATE: ln_gamma_poly f64 14 2 -0.5 0.50001
    const K1: f64 = f64::from_bits(0x3FDB0EE607209296); // 4.227843350984498e-1
    const K2: f64 = f64::from_bits(0x3FD4A34CC4A627EE); // 3.2246703342445826e-1
    const K3: f64 = f64::from_bits(0xBFB13E001A4DDB58); // -6.735230104628209e-2
    const K4: f64 = f64::from_bits(0x3F951322ABC8D183); // 2.05808083866983e-2
    const K5: f64 = f64::from_bits(0xBF7E404FE0D2647B); // -7.385551475772019e-3
    const K6: f64 = f64::from_bits(0x3F67ADD7AE79B856); // 2.8905117540523424e-3
    const K7: f64 = f64::from_bits(0xBF538ABA30F91757); // -1.1927431368276753e-3
    const K8: f64 = f64::from_bits(0x3F40B33C50A812BD); // 5.096477992202125e-4
    const K9: f64 = f64::from_bits(0xBF2D43F42126BDE9); // -2.2327761742023326e-4
    const K10: f64 = f64::from_bits(0x3F1A1DD100788FC8); // 9.962642808908275e-5
    const K11: f64 = f64::from_bits(0xBF072E4DCA353C59); // -4.42140092454738e-5
    const K12: f64 = f64::from_bits(0x3EF4C7805747B82D); // 1.9816686703970507e-5
    const K13: f64 = f64::from_bits(0xBEE80017452BB114); // -1.1444261109848831e-5
    const K14: f64 = f64::from_bits(0x3ED8224192C0E921); // 5.7539494200903685e-6

    let f2 = f * f;
    let f4 = f2 * f2;
    let f8 = f4 * f4;
    f * ((((K1 + f * K2) + f2 * (K3 + f * K4)) + f4 * ((K5 + f * K6) + f2 * (K7 + f * K8)))
        + f8 * (((K9 + f * K10) + f2 * (K11 + f * K12)) + f4 * (K13 + f * K14)))
}
