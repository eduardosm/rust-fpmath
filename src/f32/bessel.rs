//! Bessel functions for `f32`.
//!
//! They use the same formulas as the `f64` functions (see
//! `crate::f64::bessel`), but evaluated with plain `f64` arithmetic:
//!
//! * With small arguments, shorter polynomials of `x^2`.
//! * With medium arguments, the polynomials of the tables of the `f64`
//!   functions (see `ZeroPoly::eval_short`). The centers of their intervals
//!   are `f64` values that are much closer to the zeros than any `f32`, so
//!   the results close to a zero are as accurate as the others.
//! * With large arguments, the polynomials of the modulus and the phase of
//!   the `f64` functions, where the phase is still calculated as a sum of two
//!   `f64` (see `Large::eval_short`), which makes it accurate enough for the
//!   `f32` values that are closest to a zero.
//!
//! `jn` and `yn` with other orders use the recurrences of the `f64`
//! functions. For `x < n`, where the functions do not oscillate, they are
//! calculated with plain `f64` arithmetic (see `jn_backward` and
//! `yn_forward`), unless the order is large (see `ORDER_PLAIN_MAX`).
//! Otherwise, the `f64` functions are used: close to a zero, the terms of the
//! recurrence cancel each other (and the phase of the asymptotic expansions of
//! the large orders is a large number), which needs double-word arithmetic
//! also for an `f32` result. They detect the results that are zero or
//! infinite as `f32` before calculating them, with the limits of `f32` (see
//! `jn_finite` and `yn_finite` in `crate::f64::bessel`).

use super::log::ln_f64_accurate;
use crate::f64::{
    J0_LARGE_MIN, J0_SMALL_MAX, J1_LARGE_MIN, J1_SMALL_MAX, ORDER_DEBYE_MIN, Y0_LARGE_MIN,
    Y0_SMALL_MAX, Y1_LARGE_MIN, Y1_SMALL_MAX, j0_large_short, j0_zero_poly, j1_large_short,
    j1_zero_poly, jn_convergent_terms, jn_underflows, ln_1p_q, y0_large_short, y0_y1_large_short,
    y0_zero_poly, y1_large_short, y1_zero_poly, yn_overflows,
};
use crate::traits::Float as _;

// GENERATE: consts f64 FRAC_2_PI
const FRAC_2_PI: f64 = f64::from_bits(0x3FE45F306DC9C883); // 6.366197723675814e-1

/// Orders below which `jn(n, x)` and `yn(n, x)` are calculated with plain
/// `f64` arithmetic for `x < n`, which is faster than the recurrences of the
/// `f64` functions. From it, those use asymptotic expansions, which take a
/// time that does not depend on the order, for `x >= DEBYE_BELOW_MIN * n`,
/// and the results with smaller arguments are zero (`jn`) or infinite (`yn`)
/// as `f32`, which they detect beforehand (see `ORDER_DEBYE_MIN` in
/// `crate::f64::bessel`). The latter holds from the order 140 (see
/// `test_plain_max`).
///
/// The recurrences use a rounded `2 / x`, whose error is equivalent to a
/// relative change of up to 2^-53 in `x`. That changes the results by up to
/// about `n * 2^-52` relative to them, which is less than 2^-40 with these
/// orders.
const ORDER_PLAIN_MAX: u32 = ORDER_DEBYE_MIN;

// The error bounds that refer to `ORDER_PLAIN_MAX` hold with these orders,
// and so do the results that it says are zero or infinite
const _: () = assert!(140 <= ORDER_PLAIN_MAX && ORDER_PLAIN_MAX <= 1 << 12);

/// Returns `true` if `jn(n, x)` and `yn(n, x)` are calculated with plain
/// `f64` arithmetic (see `ORDER_PLAIN_MAX`), for `n >= 2` and finite `x > 0`.
#[inline]
fn uses_plain(n: u32, x: f32) -> bool {
    n < ORDER_PLAIN_MAX && x < n as f32
}

/// Orders from which `jn(n, x)` and `yn(n, x)` check whether their results
/// are zero or infinite as `f32` before the recurrences with plain `f64`
/// arithmetic (see `jn_underflows` and `yn_overflows` in
/// `crate::f64::bessel`), which would only find it out after many steps.
///
/// That check takes about as long as 30 steps of the recurrences, which is
/// less than 6% of them with these orders.
const ORDER_CHECKED_MIN: u32 = 1 << 9;

impl crate::generic::Bessel for f32 {
    #[inline]
    fn j0_finite(x: Self) -> Self {
        j0_f64(x.abs()) as f32
    }

    #[inline]
    fn j1_finite(x: Self) -> Self {
        let r = j1_f64(x.abs()) as f32;
        if x.is_sign_negative() { -r } else { r }
    }

    #[inline]
    fn y0_finite(x: Self) -> Self {
        y0_f64(x) as f32
    }

    #[inline]
    fn y1_finite(x: Self) -> Self {
        y1_f64(x) as f32
    }

    #[inline]
    fn jn_finite(n: u32, x: Self) -> Self {
        let xd = f64::from(x);
        let r = if uses_plain(n, x) {
            // The results below 2^-150 are rounded to zero as `f32`
            if n >= ORDER_CHECKED_MIN && jn_underflows(n, xd, 150) {
                0.0
            } else {
                jn_backward(n, x)
            }
        } else {
            crate::f64::jn_finite(n, xd, 150)
        };
        r as f32
    }

    #[inline]
    fn yn_finite(n: u32, x: Self) -> Self {
        let xd = f64::from(x);
        let r = if uses_plain(n, x) {
            // The results below -2^128 are rounded to an infinity as `f32`
            if n >= ORDER_CHECKED_MIN && yn_overflows(n, xd, 128) {
                f64::NEG_INFINITY
            } else {
                yn_forward(n, x)
            }
        } else {
            crate::f64::yn_finite(n, xd, 128)
        };
        r as f32
    }
}

/// Returns `j0(x)`, for finite `x >= 0`.
///
/// The comparisons use `f32` values and each branch converts its argument
/// to `f64` (see `super::erf`).
#[inline]
fn j0_f64(x: f32) -> f64 {
    if x < J0_SMALL_MAX as f32 {
        j0_small(f64::from(x))
    } else if x < J0_LARGE_MIN as f32 {
        let x = f64::from(x);
        j0_zero_poly(x).eval_short(x)
    } else {
        j0_large_short(f64::from(x))
    }
}

/// Returns `j1(x)`, for finite `x >= 0`.
#[inline]
fn j1_f64(x: f32) -> f64 {
    if x < J1_SMALL_MAX as f32 {
        j1_small(f64::from(x))
    } else if x < J1_LARGE_MIN as f32 {
        let x = f64::from(x);
        j1_zero_poly(x).eval_short(x)
    } else {
        j1_large_short(f64::from(x))
    }
}

/// Returns `y0(x)`, for finite `x > 0`.
#[inline]
fn y0_f64(x: f32) -> f64 {
    if x < Y0_SMALL_MAX as f32 {
        y0_small(f64::from(x))
    } else if x < Y0_LARGE_MIN as f32 {
        let xd = f64::from(x);
        let (poly, anchored) = y0_zero_poly(xd);
        let g = poly.eval_short(xd);
        if anchored {
            // y0(x) = g(x) + 2/π * j0(x) * ln(x / c)
            g + FRAC_2_PI * (j0_f64(x) * ln_ratio(xd, poly.center()))
        } else {
            g
        }
    } else {
        y0_large_short(f64::from(x))
    }
}

/// Returns `y1(x)`, for finite `x > 0`.
#[inline]
fn y1_f64(x: f32) -> f64 {
    if x < Y1_SMALL_MAX as f32 {
        y1_small(f64::from(x))
    } else if x < Y1_LARGE_MIN as f32 {
        let xd = f64::from(x);
        let (poly, anchored) = y1_zero_poly(xd);
        let g = poly.eval_short(xd);
        if anchored {
            // y1(x) = g(x) / x + 2/π * j1(x) * ln(x / c)
            g / xd + FRAC_2_PI * (j1_f64(x) * ln_ratio(xd, poly.center()))
        } else {
            g
        }
    } else {
        y1_large_short(f64::from(x))
    }
}

/// Returns `jn(n, x)`, for `2 <= n < ORDER_PLAIN_MAX` and `0 < x < n`, with a
/// relative error of less than about 2^-40 (see `ORDER_PLAIN_MAX`).
///
/// It uses the backward recurrence of the `f64` function (see `jn_backward`
/// in `crate::f64::bessel`), with plain `f64` arithmetic:
///
/// `a[j-1] = (2 * j / x) * a[j] - a[j+1]`
///
/// from `a[m+1] = 0` and `a[m] = 1`, where the result is `a[n]` normalized
/// with `j0(x)` or `j1(x)` (the one with the largest magnitude).
fn jn_backward(n: u32, x: f32) -> f64 {
    // m = n + k makes the relative error of the result about `c / q[k]^2`,
    // where `c` is up to about `n^(2/3) / 8 < 2^7` (see `jn_backward_terms`
    // in `crate::f64::bessel`), so it is less than 2^-49
    let k = jn_convergent_terms(n, f64::from(x), f64::exp2i_fast(28));

    // a[j + 1] and a[j], down to a[n]
    let h = 2.0 / f64::from(x);
    let (mut b, mut a) = (0.0, 1.0);
    let mut jf = f64::from(n + k);
    for _ in 0..k {
        let c = (h * jf) * a - b;
        b = a;
        a = c;
        jf -= 1.0;
    }
    let an = a;

    // a[j] / a[n] = J[j](x) / J[n](x) and |J[j](x)| <= 1, so the result is
    // less than 2^-160, which rounds to zero as `f32`, when the magnitude of
    // `a` goes beyond this limit. So `a` cannot overflow and it does not
    // need the scaling of the `f64` function: a[n] is `2 * m / x` with
    // k = 1 and less than `q[k] < 2^28 * (2 * m / x)` otherwise, which is
    // less than 2^193, and each step multiplies the magnitude by up to
    // `2 * j / x < 2^165`.
    let a_max = an * f64::exp2i_fast(160);
    for _ in 0..n {
        let c = (h * jf) * a - b;
        b = a;
        a = c;
        jf -= 1.0;
        if a.abs() > a_max {
            return 0.0;
        }
    }

    // `a` and `b` are proportional to `j0(x)` and `j1(x)`
    if a.abs() >= b.abs() {
        j0_f64(x) * (an / a)
    } else {
        j1_f64(x) * (an / b)
    }
}

/// Returns `yn(n, x)`, for `2 <= n < ORDER_PLAIN_MAX` and `0 < x < n`, with a
/// relative error of less than about 2^-40 (see `ORDER_PLAIN_MAX`), or
/// negative infinity when it is less than -2^128 (so it is rounded to an
/// infinity as `f32`).
///
/// It uses the recurrence of the `f64` function (see `yn_finite` in
/// `crate::f64::bessel`) from `y0(x)` and `y1(x)`, with plain `f64`
/// arithmetic:
///
/// `Y[k+1](x) = (2 * k / x) * Y[k](x) - Y[k-1](x)`
///
/// With these arguments, the function is negative and its magnitude is
/// close to its modulus, so the errors of the recurrence, which are small
/// relative to the modulus, are also small relative to the result.
fn yn_forward(n: u32, x: f32) -> f64 {
    let h = 2.0 / f64::from(x);
    let (mut a, mut b) = if x >= Y0_LARGE_MIN as f32 {
        // Both are calculated from their modulus and phase (`Y1_LARGE_MIN` is
        // smaller), which share some parts
        y0_y1_large_short(f64::from(x))
    } else {
        (y0_f64(x), y1_f64(x))
    };
    // With the orders `k <= x`, the magnitudes are at most about 1, so the
    // ones beyond 2^128 are reached with `k > x`, where the function is
    // negative and `|Y[k+1]| >= 2 * |Y[k]| - |Y[k-1]|`, so they keep growing
    // with the order
    let overflow = -f64::exp2i_fast(128);
    for k in 1..n {
        if b < overflow {
            return f64::NEG_INFINITY;
        }
        let c = (h * f64::from(k)) * b - a;
        a = b;
        b = c;
    }
    b
}

/// Returns `j0(x)`, for `|x| < J0_SMALL_MAX`, with a relative error of about
/// 2^-50.
#[inline]
fn j0_small(x: f64) -> f64 {
    // j0(x) ~= 1 + K1 * x^2 + ... + K6 * x^12
    // GENERATE: bessel_small_poly j0 f64 6 1
    const K1: f64 = f64::from_bits(0xBFCFFFFFFFFFFCF6); // -2.499999999999784e-1
    const K2: f64 = f64::from_bits(0x3F8FFFFFFFFC66CA); // 1.56249999995909e-2
    const K3: f64 = f64::from_bits(0xBF3C71C71996B076); // -4.3402777518042974e-4
    const K4: f64 = f64::from_bits(0x3EDC71C504B08514); // 6.781676414242534e-6
    const K5: f64 = f64::from_bits(0xBE72338F4C76F960); // -6.780551863644183e-8
    const K6: f64 = f64::from_bits(0x3DFFCAD5979A9B9C); // 4.6263916968555886e-10

    if x.exponent() < -27 {
        // j0(x) = 1 - x^2 / 4 + ..., where the terms after the first one are
        // less than 2^-56. Skipping them avoids the powers of `x` that are
        // subnormal (in `f64`) when `x` is a subnormal `f32`.
        return 1.0;
    }

    let u = x * x;
    let u2 = u * u;
    1.0 + u * (((K1 + u * K2) + u2 * (K3 + u * K4)) + (u2 * u2) * (K5 + u * K6))
}

/// Returns `j1(x)`, for `|x| < J1_SMALL_MAX`, with a relative error of about
/// 2^-50.
#[inline]
fn j1_small(x: f64) -> f64 {
    // j1(x) ~= x * (K0 + K1 * x^2 + ... + K8 * x^16)
    // GENERATE: bessel_small_poly j1 f64 9 2.375
    const K0: f64 = f64::from_bits(0x3FDFFFFFFFFFFFFE); // 4.999999999999999e-1
    const K1: f64 = f64::from_bits(0xBFAFFFFFFFFFFEAE); // -6.2499999999997655e-2
    const K2: f64 = f64::from_bits(0x3F6555555554F784); // 2.604166666656251e-3
    const K3: f64 = f64::from_bits(0xBF0C71C71C493E48); // -5.425347220422129e-5
    const K4: f64 = f64::from_bits(0x3EA6C16C0DE74700); // 6.781683870532389e-7
    const K5: f64 = f64::from_bits(0xBE3845C670BF31B6); // -5.6513955841004976e-9
    const K6: f64 = f64::from_bits(0x3DC27DFDE3B48767); // 3.363703453200258e-11
    const K7: f64 = f64::from_bits(0xBD45149D73E31110); // -1.4978615707608637e-13
    const K8: f64 = f64::from_bits(0x3CC17A7669CD3AD5); // 4.851222019480979e-16

    if x.exponent() < -27 {
        // j1(x) = x / 2 - x^3 / 16 + ..., where the terms after the first
        // one are less than 2^-57 relative to it (see `j0_small`)
        return 0.5 * x;
    }

    let u = x * x;
    let u2 = u * u;
    let u4 = u2 * u2;
    x * ((((K0 + u * K1) + u2 * (K2 + u * K3)) + u4 * ((K4 + u * K5) + u2 * (K6 + u * K7)))
        + (u4 * u4) * K8)
}

/// Returns `y0(x)`, for `0 < x < Y0_SMALL_MAX`, with a relative error of about
/// 2^-50.
#[inline]
fn y0_small(x: f64) -> f64 {
    // y0(x) ~= 2/π * ln(x) * j0(x) + K0 + K1 * x^2 + ... + K5 * x^10
    // GENERATE: bessel_small_poly y0 f64 6 0.5
    const K0: f64 = f64::from_bits(0xBFB2E4D699CBD009); // -7.380429510868693e-2
    const K1: f64 = f64::from_bits(0x3FC6BBCB41033BAE); // 1.776060168690185e-1
    const K2: f64 = f64::from_bits(0xBF9075B1BBEC32E0); // -1.60739680241474e-2
    const K3: f64 = f64::from_bits(0x3F41A6205CEE91BA); // 5.386026403937851e-4
    const K4: f64 = f64::from_bits(0xBEE3E97E2F0B7C30); // -9.494820433385021e-6
    const K5: f64 = f64::from_bits(0x3E7BA45D615F3CF9); // 1.0297466021481052e-7

    // Both terms are negative
    let t = FRAC_2_PI * (ln_f64_accurate(x, 0) * j0_small(x));
    if x.exponent() < -27 {
        // The terms after `K0` are less than 2^-56 (see `j0_small`)
        return t + K0;
    }

    let u = x * x;
    let u2 = u * u;
    t + (((K0 + u * K1) + u2 * (K2 + u * K3)) + (u2 * u2) * (K4 + u * K5))
}

/// Returns `y1(x)`, for `0 < x < Y1_SMALL_MAX`, with a relative error of about
/// 2^-46.
#[inline]
fn y1_small(x: f64) -> f64 {
    // y1(x) ~= 2/π * ln(x) * j1(x) - 2 / (π * x)
    //          + x * (K0 + K1 * x^2 + ... + K5 * x^10)
    // GENERATE: bessel_small_poly y1 f64 6 1
    const K0: f64 = f64::from_bits(0xBFC91866143C811C); // -1.9605709064581667e-1
    const K1: f64 = f64::from_bits(0x3FABD3975C63189F); // 5.4348688152047624e-2
    const K2: f64 = f64::from_bits(0xBF6835B970E687B6); // -2.9553052801966863e-3
    const K3: f64 = f64::from_bits(0x3F12C7D919C0DCD6); // 7.164251877688223e-5
    const K4: f64 = f64::from_bits(0xBEB0A666547951DA); // -9.924172718239717e-7
    const K5: f64 = f64::from_bits(0x3E42C4BC0545A009); // 8.739761061304174e-9

    // -2 / (π * x), which is the largest term
    let r = -FRAC_2_PI / x;
    if x.exponent() < -27 {
        // The other terms are less than 2^-50 relative to it (see
        // `j0_small`)
        return r;
    }

    // The three terms are negative
    let t = FRAC_2_PI * (ln_f64_accurate(x, 0) * j1_small(x));
    let u = x * x;
    let u2 = u * u;
    let p = ((K0 + u * K1) + u2 * (K2 + u * K3)) + (u2 * u2) * (K4 + u * K5);
    r + (t + x * p)
}

/// Returns `ln(x / c)`, for positive and normal `x` and `c`, with a relative
/// error of about 2^-45.
#[inline]
fn ln_ratio(x: f64, c: f64) -> f64 {
    // `x - c` is exact
    let q = (x - c) / c;
    if q.abs() < 1.0 / 256.0 {
        // ln(1 + q) ~= q - q^2 / 2 + q^3 * Q(q)
        q + (q * q) * (q * ln_1p_q(q) - 0.5)
    } else {
        // The mantissa of the quotient has 53 bits, so the absolute error of
        // the logarithm is about 2^-53, while its magnitude is at least
        // about 2^-8
        ln_f64_accurate(x / c, 0)
    }
}

#[cfg(test)]
mod tests {
    use super::{ORDER_PLAIN_MAX, jn_backward, jn_underflows, yn_forward, yn_overflows};
    use crate::f64::DEBYE_BELOW_MIN;
    use crate::traits::Float as _;

    #[test]
    fn test_plain_max() {
        // From the order 140, the results with the arguments below
        // `DEBYE_BELOW_MIN * n` are zero (`jn`) or infinite (`yn`) as `f32`,
        // and the checks of the `f64` functions find it out (see
        // `ORDER_PLAIN_MAX`). Their bounds get closer to the limits as the
        // argument increases (up to the order), and farther from them, in
        // proportion to the order, as it increases, so the largest argument
        // of each order is enough.
        let largest_arg = |n: u32| {
            let limit = DEBYE_BELOW_MIN * f64::from(n);
            let mut x = limit as f32;
            while f64::from(x) >= limit {
                x = f32::from_bits(x.to_bits() - 1);
            }
            f64::from(x)
        };
        let orders = (140..=1 << 12).chain([1 << 14, 1 << 20, (1 << 31) - 1, 1 << 31]);
        for n in orders {
            let x = largest_arg(n);
            assert!(jn_underflows(n, x, 150), "n = {n}, x = {x}");
            assert!(yn_overflows(n, x, 128), "n = {n}, x = {x}");
        }
        // It does not hold with the previous order
        assert!(!jn_underflows(139, largest_arg(139), 150));
    }

    /// Checks the relative errors of the recurrences with plain `f64`
    /// arithmetic, taking the results of the `f64` functions as exact.
    #[test]
    fn test_plain_recurrences() {
        for n in [2, 3, 8, 64, 300, 700, ORDER_PLAIN_MAX - 1] {
            // The error that grows with the order (see `ORDER_PLAIN_MAX`),
            // and the one of the functions of orders 0 and 1
            let max_err = f64::from(n) * f64::exp2i_fast(-52) + f64::exp2i_fast(-45);

            // Values spread over (0, n)
            let nf = n as f32;
            let count = (400_000 / n).max(2000);
            for i in 1..count {
                let x = (nf * (i as f32 / count as f32)).purify();
                let xd = f64::from(x);

                // The results that are not normal as `f32` have fewer
                // significant bits, or none
                let j = crate::f64::jn_finite(n, xd, 1075);
                if j >= f64::from(f32::MIN_POSITIVE) {
                    let err = ((jn_backward(n, x) - j) / j).abs();
                    assert!(err < max_err, "n = {n}, x = {x:e}, err = {err:e}");
                }
                let y = crate::f64::yn_finite(n, xd, 1024);
                if y >= -f64::from(f32::MAX) {
                    let err = ((yn_forward(n, x) - y) / y).abs();
                    assert!(err < max_err, "n = {n}, x = {x:e}, err = {err:e}");
                } else if y as f32 == f32::NEG_INFINITY {
                    assert_eq!(
                        yn_forward(n, x) as f32,
                        f32::NEG_INFINITY,
                        "n = {n}, x = {x:e}"
                    );
                }
            }
        }
    }
}
