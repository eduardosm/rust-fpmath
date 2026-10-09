//! Inverse trigonometric functions for `f32`.
//!
//! They use the same methods as the `f64` functions (see
//! `crate::f64::inv_trigonometric`), but evaluated with plain `f64`
//! arithmetic, with shorter polynomials and the table of `atan(i / 64)`
//! rounded to `f64`:
//! * `atan` and `atan2`: `atan(z) = atan(c) + atan(t)`, with
//!   `t = (z - c) / (1 + z * c)` and `c = i / 64` (see `atan_unit`).
//! * `asin` and `acos`: `asin(x) = x + x^3 * P(x^2)` for `|x| <= 1/2` (see
//!   `asin_small`), and `asin(|x|) = π/2 - 2 * asin(sqrt((1 - |x|) / 2))`
//!   otherwise (see `asin_half`).
//!
//! The relative error before the final rounding comes mostly from the
//! polynomial of `asin` (see `asin_poly`), amplified by up to 2 for
//! `|x| > 1/2`, or from `atan_unit` (for `atan` and `atan2`). The functions
//! in degrees and half-turns multiply the result by `180 / π` or `1 / π`
//! rounded to `f64` (with relative errors of 2^-54.7 and 2^-53.8), which adds
//! an error of up to about 2^-52.4 with the rounding of the product.

use crate::f64::{atan_index, atan_tbl, fast_sqrt};

// GENERATE: consts f64 FRAC_PI_2 PI FRAC_180_PI FRAC_1_PI
const FRAC_PI_2: f64 = f64::from_bits(0x3FF921FB54442D18); // 1.5707963267948966e0
const PI: f64 = f64::from_bits(0x400921FB54442D18); // 3.141592653589793e0
const FRAC_180_PI: f64 = f64::from_bits(0x404CA5DC1A63C1F8); // 5.729577951308232e1
const FRAC_1_PI: f64 = f64::from_bits(0x3FD45F306DC9C883); // 3.183098861837907e-1

impl crate::generic::InvTrigonometric for f32 {
    // GENERATE: consts f32 PI FRAC_PI_2
    const PI: f32 = f32::from_bits(0x40490FDB); // 3.1415927e0
    const FRAC_PI_2: f32 = f32::from_bits(0x3FC90FDB); // 1.5707964e0

    #[inline]
    fn asin_finite(x: Self) -> Self {
        asin_core(x) as f32
    }

    #[inline]
    fn acos_finite(x: Self) -> Self {
        acos_core(x) as f32
    }

    #[inline]
    fn atan_finite(x: Self) -> Self {
        atan_core(x) as f32
    }

    #[inline]
    fn atan2_finite(y: Self, x: Self) -> Self {
        atan2_core(y, x, 1.0) as f32
    }

    #[inline]
    fn asind_finite(x: Self) -> Self {
        (asin_core(x) * FRAC_180_PI) as f32
    }

    #[inline]
    fn acosd_finite(x: Self) -> Self {
        (acos_core(x) * FRAC_180_PI) as f32
    }

    #[inline]
    fn atand_finite(x: Self) -> Self {
        (atan_core(x) * FRAC_180_PI) as f32
    }

    #[inline]
    fn atan2d_finite(y: Self, x: Self) -> Self {
        atan2_core(y, x, FRAC_180_PI) as f32
    }

    #[inline]
    fn asinpi_finite(x: Self) -> Self {
        (asin_core(x) * FRAC_1_PI) as f32
    }

    #[inline]
    fn acospi_finite(x: Self) -> Self {
        (acos_core(x) * FRAC_1_PI) as f32
    }

    #[inline]
    fn atanpi_finite(x: Self) -> Self {
        (atan_core(x) * FRAC_1_PI) as f32
    }

    #[inline]
    fn atan2pi_finite(y: Self, x: Self) -> Self {
        atan2_core(y, x, FRAC_1_PI) as f32
    }
}

/// Returns `asin(x)`, for `|x| < 1`.
#[inline(always)]
fn asin_core(x: f32) -> f64 {
    let x = f64::from(x);
    let ax = x.abs();
    if ax <= 0.5 {
        asin_small(x)
    } else {
        // asin(|x|) = π/2 - 2 * asin(y), with y = sqrt((1 - |x|) / 2)
        (FRAC_PI_2 - 2.0 * asin_half(ax)).copysign(x)
    }
}

/// Returns `acos(x)`, for `|x| < 1`.
#[inline(always)]
fn acos_core(x: f32) -> f64 {
    let x = f64::from(x);
    let ax = x.abs();
    if ax <= 0.5 {
        // acos(x) = π/2 - asin(x)
        FRAC_PI_2 - asin_small(x)
    } else if x.is_sign_negative() {
        // acos(x) = π - 2 * asin(y), with y = sqrt((1 - |x|) / 2)
        PI - 2.0 * asin_half(ax)
    } else {
        // acos(x) = 2 * asin(y), with y = sqrt((1 - x) / 2)
        2.0 * asin_half(ax)
    }
}

/// Returns `asin(sqrt((1 - x) / 2))`, for `1/2 < x < 1`.
#[inline]
fn asin_half(x: f64) -> f64 {
    // z = (1 - x) / 2 is exact, and the error of its square root (see
    // `fast_sqrt`) is small compared to the one of the polynomial
    let z = (1.0 - x) * 0.5;
    let y = fast_sqrt(z).0;
    y + (y * z) * asin_poly(z)
}

/// Returns `asin(x)`, for `|x| <= 1/2`.
#[inline]
fn asin_small(x: f64) -> f64 {
    let u = x * x;
    x + (x * u) * asin_poly(u)
}

/// Returns `P(u)` such that `asin(x) ~= x + x^3 * P(x^2)`, for
/// `0 <= u = x^2 <= 1/4`, with a relative error of 2^-39.9 in `asin(x)`.
#[inline]
fn asin_poly(u: f64) -> f64 {
    // GENERATE: asin_poly f64 8 0.5
    const K1: f64 = f64::from_bits(0x3FC555555479BF78); // 1.6666666626724314e-1
    const K2: f64 = f64::from_bits(0x3FB333340BC1A94B); // 7.500005042098683e-2
    const K3: f64 = f64::from_bits(0x3FA6DB24EB7F2A75); // 4.464068770137867e-2
    const K4: f64 = f64::from_bits(0x3F9F2836DA5E808F); // 3.0426842764086855e-2
    const K5: f64 = f64::from_bits(0x3F9663412BCB6602); // 2.186300115395135e-2
    const K6: f64 = f64::from_bits(0x3F9527A6D5769F68); // 2.0659071711321003e-2
    const K7: f64 = f64::from_bits(0x3F5FEC145D6755A4); // 1.9483755946699622e-3
    const K8: f64 = f64::from_bits(0x3FA0DDC06A86515B); // 3.294183063840946e-2

    let u2 = u * u;
    let u4 = u2 * u2;
    ((K1 + u * K2) + u2 * (K3 + u * K4)) + u4 * ((K5 + u * K6) + u2 * (K7 + u * K8))
}

/// Returns `atan(x)`, for finite `x`.
#[inline]
fn atan_core(x: f32) -> f64 {
    let x = f64::from(x);
    let ax = x.abs();
    let r = if ax <= 1.0 {
        atan_unit(ax)
    } else {
        // atan(|x|) = π/2 - atan(1 / |x|), where `1 / |x|` is rounded to
        // `f64`
        FRAC_PI_2 - atan_unit(1.0 / ax)
    };
    r.copysign(x)
}

/// Returns `atan2(y, x) * k`, for finite and non-zero `x` and `y`.
#[inline]
fn atan2_core(y: f32, x: f32, k: f64) -> f64 {
    // atan2(y, x) = base + sign * atan(n / d), with n = min(|x|, |y|) and
    // d = max(|x|, |y|):
    // * x > 0, |y| <= |x|: atan2(|y|, x) = atan(n / d)
    // * x > 0, |y| > |x|: atan2(|y|, x) = π/2 - atan(n / d)
    // * x < 0, |y| <= |x|: atan2(|y|, x) = π - atan(n / d)
    // * x < 0, |y| > |x|: atan2(|y|, x) = π/2 + atan(n / d)
    // * y < 0: atan2(y, x) = -atan2(|y|, x)
    // `(base, sign)` is taken from `BASE_SIGN`, without branches.
    const BASE_SIGN: [(f64, f64); 8] = [
        (0.0, 1.0),
        (FRAC_PI_2, -1.0),
        (PI, -1.0),
        (FRAC_PI_2, 1.0),
        (-0.0, -1.0),
        (-FRAC_PI_2, 1.0),
        (-PI, 1.0),
        (-FRAC_PI_2, -1.0),
    ];

    let ay = y.abs();
    let ax = x.abs();
    let swap = ay > ax;
    let n = f64::from(ay.min(ax));
    let d = f64::from(ay.max(ax));
    let i = (usize::from(y.is_sign_negative()) << 2)
        | (usize::from(x.is_sign_negative()) << 1)
        | usize::from(swap);
    let (base, sign) = BASE_SIGN[i];

    // `n / d` is rounded to `f64`, and it cannot underflow (`f32` values are
    // greater than 2^-150)
    (base + sign * atan_unit(n / d)) * k
}

/// Returns `atan(z)`, for `0 <= z <= 1`, with a relative error of about
/// 2^-49.
///
/// `atan(z) = atan(c) + atan(t)`, with `t = (z - c) / (1 + z * c)` and
/// `c = i / 64` (so |t| <= 1/128), where `atan(c)` comes from a table and
/// `atan(t)` is approximated with a polynomial.
#[inline]
fn atan_unit(z: f64) -> f64 {
    // atan(t) ~= t + t^3 * (K2 + K4 * t^2)
    // GENERATE: atan_poly f64 2 1e-30 0.0079
    const K2: f64 = f64::from_bits(0xBFD5555555137350); // -3.3333333309365276e-1
    const K4: f64 = f64::from_bits(0x3FC99931798AC818); // 1.9998758730577104e-1

    let (i, c) = atan_index(z);
    // `z - c` is exact (Sterbenz lemma, or c = 0), and `1 + z * c` has a
    // relative error of at most 2^-53 (it is exact when `z` comes from an
    // `f32`)
    let t = (z - c) / (1.0 + z * c);
    let t2 = t * t;
    atan_tbl(i) + (t + (t * t2) * (K2 + t2 * K4))
}
