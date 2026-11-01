//! Bessel functions for `f64`, and the parts that `f32` shares.
//!
//! `j0` and `j1` are the Bessel functions of the first kind of orders 0 and
//! 1, and `y0` and `y1` the ones of the second kind. `j0(-x) = j0(x)` and
//! `j1(-x) = -j1(x)`, so they are calculated from `|x|`.
//!
//! # Small arguments
//!
//! They are approximated with polynomials of `u = x^2` (see `j0_small`,
//! `j1_small`, `y0_small` and `y1_small`):
//!
//! * `j0(x) = 1 + u * P(u)`
//! * `j1(x) = x * P(u)`
//! * `y0(x) = 2/π * ln(x) * j0(x) + P(u)`
//! * `y1(x) = 2/π * ln(x) * j1(x) - 2 / (π * x) + x * P(u)`
//!
//! # Medium arguments
//!
//! The functions oscillate, with a zero in each interval of length π, and
//! their results close to a zero have to be as accurate as the others
//! (relative to the result). So each function has a polynomial for each of
//! its zeros below `J0_LARGE_MIN` or `J1_LARGE_MIN` (see `ZeroPoly`), which
//! approximates the function in an interval of length π around the zero:
//!
//! `f(c + d) ~= K0 + K1 * d + ... + K20 * d^20`
//!
//! where `c` is a value very close to the zero, such that `d = x - c` is
//! exact in the whole interval, and `K0 = f(c)` is tiny. Hence, the smallest
//! results, `K0 + K1 * d`, do not lose accuracy.
//!
//! `y0` and `y1` have a singularity at zero, which makes the polynomials for
//! their first zeros inaccurate, so those approximate the functions without
//! the singularity instead, which are entire (see `y0_anchored` and
//! `y1_anchored`):
//!
//! * `y0(x) = g(x) + 2/π * j0(x) * ln(x / c)`
//! * `y1(x) = g(x) / x + 2/π * j1(x) * ln(x / c)`
//!
//! The logarithm of `x / c` is zero at `c`, so its term is also small close
//! to the zero.
//!
//! # Large arguments
//!
//! They are calculated from the modulus and the phase of the functions of
//! each order (see `Large`):
//!
//! * `j(x) = M(x) * cos(θ(x))` and `y(x) = M(x) * sin(θ(x))`
//! * `M(x) = sqrt(2 / (π * x)) * A(1 / x^2)`
//! * `θ(x) = x - (2 * order + 1) * π/4 + B(1 / x^2) / x`
//!
//! where `A` and `B` are approximated with polynomials, as in their
//! asymptotic series. `x` is reduced like in the trigonometric functions (see
//! `super::trigonometric`), the last term of the phase is added to the
//! reduced argument, and the sine is evaluated like in `sin`.
//!
//! The results close to a zero are about `M(x) * r`, where `r` is the distance
//! from the phase to a multiple of π (or to an odd multiple of π/2), which
//! can be as small as about 2^-60 while the terms of the phase that cancel
//! are much larger. So when `|r| < NEAR_ZERO`, the phase is calculated again
//! with integer arithmetic, with an absolute error of about 2^-130 (see
//! `phase_precise`).
//!
//! # Accuracy
//!
//! The polynomials, the logarithms, the sines and their products are
//! calculated as sums of two `f64`. The relative error before the final
//! rounding is about 2^-68 with small and medium arguments, and about 2^-65
//! with large ones (the one of the sine, see `large_eval`).
//!
//! # Other orders
//!
//! `jn` and `yn` use the recurrence
//!
//! `C[k+1](x) = (2 * k / x) * C[k](x) - C[k-1](x)`
//!
//! with double-word arithmetic. It is stable from the orders 0 and 1 for the
//! functions of the second kind, and for the ones of the first kind up to
//! the order `x` (see `jn_forward` and `yn_finite`). With larger orders, the
//! functions of the first kind use it backwards (see `jn_backward`).
//!
//! Those recurrences take a number of steps proportional to the order. As
//! `x` decreases below the order, the functions of the first kind end up
//! underflowing to zero and the ones of the second kind overflowing, which a
//! bound detects before the recurrence (see `jn_underflows` and
//! `yn_overflows`).
//!
//! When `x` is much larger than the square of the order, the functions are
//! calculated from their modulus and phase instead, like with large
//! arguments, with the first terms of their asymptotic series (see
//! `order_large_eval`).
//!
//! With large orders, they are calculated from Debye's asymptotic expansions
//! instead (see `ORDER_DEBYE_MIN`, and `DEBYE_BELOW_MIN` for `x < n`), whose
//! variables are `R = sqrt(|n^2 - x^2|)` and its ratio to the order:
//!
//! * `x < n`: `jn(n, x)` and `yn(n, x)` are exponentials of
//!   `n * (atanh(w) - w)` (with opposite signs), with `w = R / n`, times
//!   factors that depend on `R` (see `debye_below`).
//! * `x > n`: they are calculated from their modulus and phase, where the
//!   phase is `n * (v - atan(v)) - π/4`, with `v = R / n`, plus small terms
//!   (see `debye_above`).
//!
//! The one for `x > n` is also used with smaller orders, far enough beyond
//! them, where it is faster than the recurrence (see `DEBYE_FAR_R_MIN`).
//!
//! The order is up to 2^31, so `atanh` and `atan` are calculated with an
//! error of about 2^-100 (see `atanh_precise` and `atan_precise`). The
//! expansions are not valid close to `x = n`, where Olver's uniform
//! expansions are used instead, which are expressed with the Airy functions
//! of an argument `t` that is between -25.4 and 30.9 there (see `olver`).
//! They are calculated from their Taylor series at the closest point of a
//! table (see `airy`).
//!
//! The functions do not oscillate when `x` is less than the order, and the
//! relative error before the final rounding is about 2^-66. Otherwise, the
//! error is only small relative to the modulus of the functions of that
//! order (about 2^-65 times it), so the results close to a zero are less
//! accurate: the terms of the recurrence nearly cancel there, and the phase
//! is not calculated again like with orders 0 and 1.

use super::exp::{Reduced as ExpReduced, exp_poly};
use super::f64x2::{F64x2, Split};
use super::log::{ln_1p_poly_accurate, ln_accurate_parts, ln_accurate_sum_parts};
use super::trigonometric::{Reduced, mul_128_pi, reduce_rad, reduced_from_fixed};
use super::{div_residual, fast_sqrt, recip_parts, sqrt_parts, square_parts};
use crate::traits::Float as _;

// GENERATE: consts F64x2 FRAC_2_PI
const FRAC_2_PI: F64x2 = F64x2::from_bits(0x3FE45F306DC9C883, 0xBC86B01EC5417056); // 6.366197723675813430755350534901e-1

// GENERATE: consts f64 FRAC_1_PI
const FRAC_1_PI: f64 = f64::from_bits(0x3FD45F306DC9C883); // 3.183098861837907e-1

// GENERATE: split_const FRAC_PI_128_R PI -7 53 53
const FRAC_PI_128_R0: f64 = f64::from_bits(0x3F9921FB54442D18); // 2.454369260617026e-2
const FRAC_PI_128_R1: f64 = f64::from_bits(0x3C31A62633145C07); // 9.567553118338697e-19

/// `j0` is calculated as `1 + x^2 * P(x^2)` for `|x| < J0_SMALL_MAX`.
pub(crate) const J0_SMALL_MAX: f64 = 1.0;

/// `j1` is calculated as `x * P(x^2)` for `|x| < J1_SMALL_MAX`.
pub(crate) const J1_SMALL_MAX: f64 = 2.375;

/// `y0` is calculated as `2/π * ln(x) * j0(x) + P(x^2)` for
/// `x < Y0_SMALL_MAX`.
pub(crate) const Y0_SMALL_MAX: f64 = 0.5;

/// `y1` is calculated as `2/π * ln(x) * j1(x) - 2 / (π * x) + x * P(x^2)` for
/// `x < Y1_SMALL_MAX`.
pub(crate) const Y1_SMALL_MAX: f64 = 1.0;

/// `j0` is calculated from its modulus and phase for `|x| >= J0_LARGE_MIN`.
///
/// The zeros of `j0` and `y1` are close to `(k - 1/4) * π`, and the last
/// interval of their tables ends at `15.25 * π ~= 47.909`.
pub(crate) const J0_LARGE_MIN: f64 = 47.875;

/// `j1` is calculated from its modulus and phase for `|x| >= J1_LARGE_MIN`.
///
/// The zeros of `j1` and `y0` are close to `(k + 1/4) * π`, and the last
/// interval of their tables ends at `15.75 * π ~= 49.48`.
pub(crate) const J1_LARGE_MIN: f64 = 49.375;

/// `y0` is calculated from its modulus and phase for `x >= Y0_LARGE_MIN`.
pub(crate) const Y0_LARGE_MIN: f64 = J1_LARGE_MIN;

/// `y1` is calculated from its modulus and phase for `x >= Y1_LARGE_MIN`.
pub(crate) const Y1_LARGE_MIN: f64 = J0_LARGE_MIN;

/// Number of intervals of `Y0_TBL` whose polynomials approximate `y0` without
/// its singularity.
const Y0_ANCHORED: usize = 3;

/// Number of intervals of `Y1_TBL` whose polynomials approximate `y1` without
/// its singularity.
const Y1_ANCHORED: usize = 2;

/// With large arguments, the phase is calculated with `phase_precise` when
/// it is closer than this to a zero of the function.
///
/// The absolute error of the phase calculated by `Large::eval` is about
/// 2^-81, so its error relative to the distance to the zero is at most about
/// 2^-71.
const NEAR_ZERO: f64 = 1.0 / 1024.0;

/// With large arguments, the polynomials of the modulus and the phase are
/// only evaluated when the exponent of `x` is less than this: the terms
/// after the first one of each are less than 2^-70 relative to it otherwise.
const LARGE_POLY_MAX_EXP: i16 = 35;

/// With large arguments, the last term of the phase, which is less than
/// `1 / x`, is ignored from this exponent of `x`: it is much less than the
/// smallest distance from the phase to a zero of the function.
const LARGE_PHASE_MAX_EXP: i16 = 140;

/// With large arguments, `1 / x` is scaled from this exponent of `x`, so it
/// and its products are not subnormal (see `LargeArg`).
const LARGE_SCALED_MIN_EXP: i16 = 500;

// `Large::eval` applies the scaling of `LargeArg` only to the modulus from
// `LARGE_POLY_MAX_EXP`, and it skips the other terms that use `1 / x` from
// `LARGE_PHASE_MAX_EXP`, so the scaling cannot start below either of them
const _: () = assert!(
    LARGE_POLY_MAX_EXP <= LARGE_SCALED_MIN_EXP && LARGE_PHASE_MAX_EXP <= LARGE_SCALED_MIN_EXP
);

impl crate::generic::Bessel for f64 {
    #[inline]
    fn j0_finite(x: Self) -> Self {
        let (hi, lo) = j0_parts(x.abs());
        hi + lo
    }

    #[inline]
    fn j1_finite(x: Self) -> Self {
        let (hi, lo) = j1_parts(x.abs());
        let r = hi + lo;
        if x.is_sign_negative() { -r } else { r }
    }

    #[inline]
    fn y0_finite(x: Self) -> Self {
        let (hi, lo) = y0_parts(x);
        hi + lo
    }

    #[inline]
    fn y1_finite(x: Self) -> Self {
        let (hi, lo) = y1_parts(x);
        hi + lo
    }

    #[inline]
    fn jn_finite(n: u32, x: Self) -> Self {
        // The results below 2^-1075 are rounded to zero
        jn_finite(n, x, 1075)
    }

    #[inline]
    fn yn_finite(n: u32, x: Self) -> Self {
        // The results below -2^1024 are rounded to an infinity
        yn_finite(n, x, 1024)
    }
}

/// Approximation of a function around one of its zeros:
///
/// `f(c + d) ~= K0 + K1 * d + ... + K20 * d^20`
///
/// where `c` is very close to the zero (so `K0 = f(c)` is tiny) and
/// `d = x - c` is exact for every `x` of the interval of the zero.
pub(crate) struct ZeroPoly {
    c: f64,
    /// `K0, ..., K10`
    k_hi: [F64x2; 11],
    /// `K11, ..., K20`
    k_lo: [f64; 10],
}

impl ZeroPoly {
    /// Returns `(hi, lo)` such that `hi + lo ~= f(x)`, for `x` in the
    /// interval of the polynomial, with a relative error of about 2^-68:
    /// * Polynomial approximation: less than 2^-70 (2^-71 without the
    ///   rounding of the coefficients).
    /// * Evaluation of the terms of degree 11 and higher (less than 2^-16
    ///   relative to the result): ~2^-68.
    /// * Evaluation of the other terms: ~2^-72.
    ///
    /// The result is `K0` when `x` is `c`. Otherwise, `x` is not much closer
    /// to the zero than to `c`, so the sum of `K0` and the other terms does
    /// not lose accuracy.
    #[inline]
    fn eval(&self, x: f64) -> (f64, f64) {
        // exact
        let d = x - self.c;
        let [k0, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10] = self.k_hi;
        let [k11, k12, k13, k14, k15, k16, k17, k18, k19, k20] = self.k_lo;

        // The terms of degree 11 and higher are evaluated with plain `f64`
        // arithmetic
        let d2 = d * d;
        let d4 = d2 * d2;
        let v = (((k11 + d * k12) + d2 * (k13 + d * k14))
            + d4 * ((k15 + d * k16) + d2 * (k17 + d * k18)))
            + (d4 * d4) * (k19 + d * k20);

        // The next steps use exact products. The terms of even and odd
        // degrees have very different magnitudes when `d` is small, and the
        // partial sums can nearly cancel, so the sums are also exact.
        let ds = Split::new(d);
        let (h, l) = mul_add(ds, v, 0.0, k10);
        let (h, l) = mul_add(ds, h, l, k9);
        let (h, l) = mul_add(ds, h, l, k8);
        let (h, l) = mul_add(ds, h, l, k7);
        let (h, l) = mul_add(ds, h, l, k6);
        let (h, l) = mul_add(ds, h, l, k5);
        let (h, l) = mul_add(ds, h, l, k4);
        let (h, l) = mul_add(ds, h, l, k3);
        let (h, l) = mul_add(ds, h, l, k2);
        let (h, l) = mul_add(ds, h, l, k1);
        mul_add(ds, h, l, k0)
    }

    /// Returns `c`.
    #[inline]
    pub(crate) fn center(&self) -> f64 {
        self.c
    }

    /// Returns `f(x)`, for `x` in the interval of the polynomial, with a
    /// relative error of about 2^-50.
    ///
    /// It is meant for the `f32` functions: it evaluates the whole polynomial
    /// with plain `f64` arithmetic.
    #[inline]
    pub(crate) fn eval_short(&self, x: f64) -> f64 {
        // exact
        let d = x - self.c;
        let [k0, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10] = self.k_hi.map(F64x2::hi);
        let [k11, k12, k13, k14, k15, k16, k17, k18, k19, k20] = self.k_lo;

        let d2 = d * d;
        let d4 = d2 * d2;
        let d8 = d4 * d4;
        let v0 = ((k1 + d * k2) + d2 * (k3 + d * k4)) + d4 * ((k5 + d * k6) + d2 * (k7 + d * k8));
        let v1 =
            ((k9 + d * k10) + d2 * (k11 + d * k12)) + d4 * ((k13 + d * k14) + d2 * (k15 + d * k16));
        let v2 = (k17 + d * k18) + d2 * (k19 + d * k20);
        k0 + d * (v0 + d8 * (v1 + d8 * v2))
    }
}

/// Returns the polynomial of `J0_TBL` for the interval of `x`, for
/// `J0_SMALL_MAX <= x < J0_LARGE_MIN`.
#[inline]
pub(crate) fn j0_zero_poly(x: f64) -> &'static ZeroPoly {
    // The interval of the k-th zero is `[(k - 3/4) * π, (k + 1/4) * π)`
    &J0_TBL[(x * FRAC_1_PI + 0.75) as usize - 1]
}

/// Returns the polynomial of `J1_TBL` for the interval of `x`, for
/// `J1_SMALL_MAX <= x < J1_LARGE_MIN`.
#[inline]
pub(crate) fn j1_zero_poly(x: f64) -> &'static ZeroPoly {
    // The interval of the k-th zero is `[(k - 1/4) * π, (k + 3/4) * π)`
    &J1_TBL[(x * FRAC_1_PI + 0.25) as usize - 1]
}

/// Returns the polynomial of `Y0_TBL` for the interval of `x`, and whether
/// it approximates `y0` without its singularity (see `y0_anchored`), for
/// `Y0_SMALL_MAX <= x < Y0_LARGE_MIN`.
#[inline]
pub(crate) fn y0_zero_poly(x: f64) -> (&'static ZeroPoly, bool) {
    // The interval of the k-th zero is `[(k - 5/4) * π, (k - 1/4) * π)`
    let i = (x * FRAC_1_PI + 1.25) as usize - 1;
    (&Y0_TBL[i], i < Y0_ANCHORED)
}

/// Returns the polynomial of `Y1_TBL` for the interval of `x`, and whether
/// it approximates `y1` without its singularity (see `y1_anchored`), for
/// `Y1_SMALL_MAX <= x < Y1_LARGE_MIN`.
#[inline]
pub(crate) fn y1_zero_poly(x: f64) -> (&'static ZeroPoly, bool) {
    // The interval of the k-th zero is `[(k - 3/4) * π, (k + 1/4) * π)`
    let i = (x * FRAC_1_PI + 0.75) as usize - 1;
    (&Y1_TBL[i], i < Y1_ANCHORED)
}

/// Returns `(hi, lo)` such that `hi + lo ~= c + f * (a_hi + a_lo)`, with the
/// product calculated by `Split::mul`.
///
/// Unlike `Split::mul_add`, it does not require the product to be smaller
/// than `c`.
#[inline]
fn mul_add(f: Split, a_hi: f64, a_lo: f64, c: F64x2) -> (f64, f64) {
    let (p_hi, p_lo) = f.mul(a_hi, a_lo);
    let s = F64x2::add11(c.hi(), p_hi);
    (s.hi(), s.lo() + (p_lo + c.lo()))
}

/// Returns `(hi, lo)` such that `hi + lo ~= (a_hi + a_lo) * (b_hi + b_lo)`,
/// with a relative error of about
/// `2^-76 + 2^-52 * (|a_lo / a_hi| + |b_lo / b_hi|)`.
#[inline]
fn mul22(a_hi: f64, a_lo: f64, b_hi: f64, b_lo: f64) -> (f64, f64) {
    let (hi, lo) = Split::new(a_hi).mul(b_hi, b_lo);
    (hi, lo + a_lo * (b_hi + b_lo))
}

/// Returns `(hi, lo)` such that `hi + lo ~= (a_hi + a_lo) * 2/π`, with the
/// relative error of `mul22`.
#[inline]
fn mul_frac_2_pi(a_hi: f64, a_lo: f64) -> (f64, f64) {
    mul22(a_hi, a_lo, FRAC_2_PI.hi(), FRAC_2_PI.lo())
}

/// Returns `(hi, lo)` such that `hi + lo ~= j0(x)`, for finite `x >= 0`.
#[inline]
fn j0_parts(x: f64) -> (f64, f64) {
    if x < J0_SMALL_MAX {
        j0_small(x)
    } else if x < J0_LARGE_MIN {
        j0_zero_poly(x).eval(x)
    } else {
        large_eval(x, &LARGE0, 64)
    }
}

/// Returns `(hi, lo)` such that `hi + lo ~= j1(x)`, for finite `x >= 0`.
#[inline]
fn j1_parts(x: f64) -> (f64, f64) {
    if x.exponent() < -40 {
        // j1(x) = x / 2 - x^3 / 16 + ..., where the terms after the first one
        // are less than 2^-83 relative to it. The product is exact unless it
        // is subnormal, where it is rounded only once.
        (0.5 * x, 0.0)
    } else if x < J1_SMALL_MAX {
        j1_small(x)
    } else if x < J1_LARGE_MIN {
        j1_zero_poly(x).eval(x)
    } else {
        large_eval(x, &LARGE1, 64)
    }
}

/// Returns `(hi, lo)` such that `hi + lo ~= y0(x)`, for finite `x > 0`.
#[inline]
fn y0_parts(x: f64) -> (f64, f64) {
    if x < Y0_SMALL_MAX {
        y0_small(x)
    } else if x < Y0_LARGE_MIN {
        let (poly, anchored) = y0_zero_poly(x);
        if anchored {
            y0_anchored(x, poly)
        } else {
            poly.eval(x)
        }
    } else {
        large_eval(x, &LARGE0, 0)
    }
}

/// Returns `(hi, lo)` such that `hi + lo ~= y1(x)`, for finite `x > 0`.
#[inline]
fn y1_parts(x: f64) -> (f64, f64) {
    if x.exponent() < Y1_TINY_EXP {
        y1_tiny(x)
    } else if x < Y1_SMALL_MAX {
        y1_small(x)
    } else if x < Y1_LARGE_MIN {
        let (poly, anchored) = y1_zero_poly(x);
        if anchored {
            y1_anchored(x, poly)
        } else {
            poly.eval(x)
        }
    } else {
        large_eval(x, &LARGE1, 0)
    }
}

/// Returns `(hi, lo)` such that `hi + lo ~= j0(x)`, for `|x| < J0_SMALL_MAX`,
/// with a relative error of about 2^-71:
/// * Polynomial approximation: 2^-82.
/// * Evaluation of the terms of degree 4 and higher of the polynomial of
///   `x^2` (less than 2^-17 relative to the result): ~2^-69.
/// * Evaluation of the other terms: ~2^-74.
#[inline]
fn j0_small(x: f64) -> (f64, f64) {
    // j0(x) ~= 1 + K1 * x^2 + ... + K9 * x^18
    // GENERATE: bessel_small_poly j0 F64x2:3,f64 9 1
    const K1: F64x2 = F64x2::from_bits(0xBFD0000000000000, 0x3B38D260E1ADF8EA); // -2.499999999999999999999794678959e-1
    const K2: F64x2 = F64x2::from_bits(0x3F90000000000000, 0xBB8E7F091A03D7CA); // 1.562499999999999999919277144361e-2
    const K3: F64x2 = F64x2::from_bits(0xBF3C71C71C71C71C, 0xBBCEFF67755B59A2); // -4.340277777777777668124210726075e-4
    const K4: f64 = f64::from_bits(0x3EDC71C71C71C6C6); // 6.7816840277777046e-6
    const K5: f64 = f64::from_bits(0xBE723456789A6B17); // -6.781684027750069e-8
    const K6: f64 = f64::from_bits(0x3E002E85C02BEE25); // 4.70950279072434e-10
    const K7: f64 = f64::from_bits(0xBD8522A3BB2E7582); // -2.402806653606472e-12
    const K8: f64 = f64::from_bits(0x3D05223395ADF55F); // 9.385203556461052e-15
    const K9: f64 = f64::from_bits(0xBC807DEF0CCF5444); // -2.860893616471852e-17

    if x.exponent() < -40 {
        // j0(x) = 1 - x^2 / 4 + ..., where the terms after the first one are
        // less than 2^-82. Skipping them avoids the powers of `x` that are
        // subnormal, which can be very slow on some CPUs.
        return (1.0, 0.0);
    }

    // x^2 = u + u_lo, where `u` is `x^2` rounded
    let (u, u_lo) = square_parts(x);

    // The terms of degree 4 and higher are evaluated with plain `f64`
    // arithmetic
    let u2 = u * u;
    let w = ((K4 + u * K5) + u2 * (K6 + u * K7)) + (u2 * u2) * (K8 + u * K9);

    // The next steps use exact products
    let us = Split::new(u);
    let (h, l) = us.mul_add(w, 0.0, K3);
    let (h, l) = us.mul_add(h, l, K2);
    let (h, l) = us.mul_add(h, l, K1);
    let (p_hi, p_lo) = us.mul_add(h, l, F64x2::new1(1.0));

    // With Q(u) = u * P(u), Q(u + u_lo) ~= Q(u) + u_lo * Q'(u), where
    // Q'(u) ~= K1 + 2 * K2 * u + 3 * K3 * u^2 (the next term is less than
    // 2^-15 relative to `Q'(u)`, and |u_lo * Q'(u)| < 2^-55)
    let dq = K1.hi() + u * (2.0 * K2.hi() + u * (3.0 * K3.hi()));
    (p_hi, p_lo + u_lo * dq)
}

/// Returns `(hi, lo)` such that `hi + lo ~= j1(x)`, for
/// `2^-40 <= |x| < J1_SMALL_MAX`, with a relative error of about 2^-70:
/// * Polynomial approximation: 2^-78.
/// * Evaluation of the terms of degree 6 and higher of the polynomial of
///   `x^2` (less than 2^-17.5 relative to the result): ~2^-70.
/// * Evaluation of the other terms: ~2^-73.
#[inline]
fn j1_small(x: f64) -> (f64, f64) {
    // j1(x) ~= x * (K0 + K1 * x^2 + ... + K11 * x^22)
    // GENERATE: bessel_small_poly j1 F64x2:6,f64 12 2.375
    const K0: F64x2 = F64x2::from_bits(0x3FE0000000000000, 0xBB00CC7F917385F1); // 4.999999999999999999999982630426e-1
    const K1: F64x2 = F64x2::from_bits(0xBFB0000000000000, 0x3B593945717D69B4); // -6.249999999999999999991654172502e-2
    const K2: F64x2 = F64x2::from_bits(0x3F65555555555555, 0x3C053C06EDB46A24); // 2.604166666666666665996815310398e-3
    const K3: F64x2 = F64x2::from_bits(0xBF0C71C71C71C71C, 0xBB90D3F1B3D7F1A2); // -5.425347222222222010141149117583e-5
    const K4: F64x2 = F64x2::from_bits(0x3EA6C16C16C16BF6, 0xBB3965C4F79BB1D5); // 6.781684027777742886400826403069e-7
    const K5: F64x2 = F64x2::from_bits(0xBE3845C8A0CE4117, 0xBAD64A5564379B17); // -5.651403356478078467908977748263e-9
    const K6: f64 = f64::from_bits(0x3DC27E4FB773A270); // 3.363930569122966e-11
    const K7: f64 = f64::from_bits(0xBD4522A43D59A5D9); // -1.50175470978996e-13
    const K8: f64 = f64::from_bits(0x3CC2C974FDFB2B90); // 5.214423728357589e-16
    const K9: f64 = f64::from_bits(0xBC3AB7EA63C9EBE9); // -1.4484084642662379e-18
    const K10: f64 = f64::from_bits(0x3BAF0B4FEE5EB03B); // 3.2869312996884458e-21
    const K11: f64 = f64::from_bits(0xBB1C87B27B942C7E); // -5.899879611599825e-24

    // x^2 = u + u_lo, where `u` is `x^2` rounded
    let (u, u_lo) = square_parts(x);

    // The terms of degree 6 and higher are evaluated with plain `f64`
    // arithmetic
    let u2 = u * u;
    let w = ((K6 + u * K7) + u2 * (K8 + u * K9)) + (u2 * u2) * (K10 + u * K11);

    // The next steps use exact products
    let us = Split::new(u);
    let (h, l) = us.mul_add(w, 0.0, K5);
    let (h, l) = us.mul_add(h, l, K4);
    let (h, l) = us.mul_add(h, l, K3);
    let (h, l) = us.mul_add(h, l, K2);
    let (h, l) = us.mul_add(h, l, K1);
    let (p_hi, p_lo) = us.mul_add(h, l, K0);

    // P(u + u_lo) ~= P(u) + u_lo * P'(u), where
    // P'(u) ~= K1 + 2 * K2 * u + ... + 6 * K6 * u^5 (the next terms are less
    // than 2^-16 relative to `P'(u)`, and |u_lo * P'(u)| < 2^-53 * P(u))
    let dp = (K1.hi() + u * (2.0 * K2.hi()))
        + u2 * ((3.0 * K3.hi() + u * (4.0 * K4.hi())) + u2 * (5.0 * K5.hi() + u * (6.0 * K6)));
    let p_lo = p_lo + u_lo * dp;

    // exact products
    Split::new(x).mul(p_hi, p_lo)
}

/// Returns `(hi, lo)` such that `hi + lo ~= y0(x)`, for
/// `0 < x < Y0_SMALL_MAX`, with a relative error of about 2^-70:
/// * Polynomial approximation: 2^-76.
/// * Logarithm: 2^-75 (see `ln_accurate_parts`).
/// * `j0(x)`: 2^-71 (see `j0_small`).
/// * Evaluation of the polynomial and the products: ~2^-72.
///
/// The two terms of the sum are negative.
#[inline]
fn y0_small(x: f64) -> (f64, f64) {
    // y0(x) ~= 2/π * ln(x) * j0(x) + K0 + K1 * x^2 + ... + K7 * x^14
    // GENERATE: bessel_small_poly y0 F64x2:4,f64 8 0.5
    const K0: F64x2 = F64x2::from_bits(0xBFB2E4D699CBD01F, 0x3C5DDFE0C8690775); // -7.380429510868722527431550966996e-2
    const K1: F64x2 = F64x2::from_bits(0x3FC6BBCB41034286, 0x3C6D8F011B518D00); // 1.776060168690671420791774349967e-1
    const K2: F64x2 = F64x2::from_bits(0xBF9075B1BBF41364, 0x3C310DD98136DCE4); // -1.607396802593842507187490222271e-2
    const K3: F64x2 = F64x2::from_bits(0x3F41A6206B7B96B0, 0x3BD301D484C93096); // 5.386026668616397139480776843801e-4
    const K4: f64 = f64::from_bits(0xBEE3E997941E4A20); // -9.495005205005741e-6
    const K5: f64 = f64::from_bits(0x3E7BCE4A58646313); // 1.0358475863550803e-7
    const K6: f64 = f64::from_bits(0xBE0A6ED69771AC8D); // -7.69300441828892e-10
    const K7: f64 = f64::from_bits(0x3D92255E00D9DA1A); // 4.125914899053552e-12

    // t_hi + t_lo = 2/π * ln(x) * j0(x)
    let (xn, edelta) = x.normalize_arg();
    let (l_hi, l_lo) = ln_accurate_parts(xn, edelta);
    let (j_hi, j_lo) = j0_small(x);
    let (a_hi, a_lo) = mul22(l_hi, l_lo, j_hi, j_lo);
    let (t_hi, t_lo) = mul_frac_2_pi(a_hi, a_lo);

    if x.exponent() < -40 {
        // The terms after `K0` are less than 2^-80, and skipping them avoids
        // the powers of `x` that are subnormal
        let s = F64x2::fast_add11(t_hi, K0.hi());
        return (s.hi(), s.lo() + (t_lo + K0.lo()));
    }

    // x^2 = u + u_lo, where `u` is `x^2` rounded
    let (u, u_lo) = square_parts(x);

    // The terms of degree 4 and higher (less than 2^-21 relative to the
    // result) are evaluated with plain `f64` arithmetic
    let w = (K4 + u * K5) + (u * u) * (K6 + u * K7);

    // The next steps use exact products
    let us = Split::new(u);
    let (h, l) = us.mul_add(w, 0.0, K3);
    let (h, l) = us.mul_add(h, l, K2);
    let (h, l) = us.mul_add(h, l, K1);
    let (p_hi, p_lo) = us.mul_add(h, l, K0);
    // P(u + u_lo) ~= P(u) + u_lo * P'(u), where P'(u) ~= K1 + 2 * K2 * u
    // (the next terms are less than 2^-11 relative to `P'(u)`, and
    // |u_lo * P'(u)| is less than 2^-56 relative to the result)
    let p_lo = p_lo + u_lo * (K1.hi() + u * (2.0 * K2.hi()));

    let s = F64x2::add11(t_hi, p_hi);
    (s.hi(), s.lo() + (t_lo + p_lo))
}

/// `y1(x)` is calculated as `-2 / (π * x)` when the exponent of `x` is less
/// than this.
const Y1_TINY_EXP: i16 = -40;

/// Returns `(hi, lo)` such that `hi + lo ~= y1(x)`, for
/// `0 < x < 2^Y1_TINY_EXP`.
#[inline]
fn y1_tiny(x: f64) -> (f64, f64) {
    // y1(x) = -2 / (π * x) + x / π * ln(x) + ..., where the terms after the
    // first one are less than 2^-74 relative to it
    if x >= f64::MIN_POSITIVE {
        // The result is at most about 2^1021.4
        let r = -(FRAC_2_PI * F64x2::div11(1.0, x));
        (r.hi(), r.lo())
    } else {
        // x = m * 2^e, with 1 <= m < 2, so the reciprocal of `m` does not
        // overflow, and the result (which overflows when x < ~2^-1024.65) is
        // rounded only once
        let (xn, edelta) = x.normalize_arg();
        let e = xn.exponent() + edelta;
        let (t_hi, t_lo) = recip_parts(xn.set_exp(0));
        let (hi, lo) = mul_frac_2_pi(t_hi, t_lo);
        (-F64x2::fast_add11(hi, lo).scalbn_to_f64(-i32::from(e)), 0.0)
    }
}

/// Returns `(hi, lo)` such that `hi + lo ~= y1(x)`, for
/// `2^Y1_TINY_EXP <= x < Y1_SMALL_MAX`, with a relative error of about
/// 2^-70:
/// * Polynomial approximation: 2^-77.
/// * Reciprocal: 2^-77 (see `recip_parts`).
/// * Logarithm: 2^-75 (see `ln_accurate_parts`), in a term that is less than
///   a sixth of the result.
/// * `j1(x)`: 2^-70 (see `j1_small`), in the same term.
/// * Evaluation of the polynomial and the products: ~2^-72.
///
/// The three terms of the sum are negative.
#[inline]
fn y1_small(x: f64) -> (f64, f64) {
    // y1(x) ~= 2/π * ln(x) * j1(x) - 2 / (π * x)
    //          + x * (K0 + K1 * x^2 + ... + K8 * x^16)
    // GENERATE: bessel_small_poly y1 F64x2:4,f64 9 1
    const K0: F64x2 = F64x2::from_bits(0xBFC91866143CBC8A, 0xBC61D78B53B38D0A); // -1.960570906462389484054214186310e-1
    const K1: F64x2 = F64x2::from_bits(0x3FABD3975C75B4A7, 0x3C2EA87980AA75D7); // 5.434868816051024398127030275075e-2
    const K2: F64x2 = F64x2::from_bits(0xBF6835B97894BE5A, 0x3BD0758A999D4E16); // -2.955305336079833304566331154734e-3
    const K3: F64x2 = F64x2::from_bits(0x3F12C7DBFFCDE8C2, 0x3BA04CCA095E9032); // 7.164268749973708962758013825659e-5
    const K4: f64 = f64::from_bits(0xBEB0A780AC76B9C7); // -9.926740619326759e-7
    const K5: f64 = f64::from_bits(0x3E432E5A4D0BC1A4); // 8.931879598326294e-9
    const K6: f64 = f64::from_bits(0xBDCF0CE21DE93ECF); // -5.648021264127144e-11
    const K7: f64 = f64::from_bits(0x3D52A45EA17DCE20); // 2.64919733464964e-13
    const K8: f64 = f64::from_bits(0xBCD108B20EF3BEE7); // -9.45575115744425e-16

    // r_hi + r_lo = 2 / (π * x)
    let (t_hi, t_lo) = recip_parts(x);
    let (r_hi, r_lo) = mul_frac_2_pi(t_hi, t_lo);

    // t_hi + t_lo = 2/π * ln(x) * j1(x)
    let (l_hi, l_lo) = ln_accurate_parts(x, 0);
    let (j_hi, j_lo) = j1_small(x);
    let (a_hi, a_lo) = mul22(l_hi, l_lo, j_hi, j_lo);
    let (t_hi, t_lo) = mul_frac_2_pi(a_hi, a_lo);

    // x^2 = u + u_lo, where `u` is `x^2` rounded
    let (u, u_lo) = square_parts(x);

    // The terms of degree 4 and higher (less than 2^-19 relative to the
    // result) are evaluated with plain `f64` arithmetic
    let u2 = u * u;
    let w = ((K4 + u * K5) + u2 * (K6 + u * K7)) + (u2 * u2) * K8;

    // The next steps use exact products
    let us = Split::new(u);
    let (h, l) = us.mul_add(w, 0.0, K3);
    let (h, l) = us.mul_add(h, l, K2);
    let (h, l) = us.mul_add(h, l, K1);
    let (p_hi, p_lo) = us.mul_add(h, l, K0);
    // P(u + u_lo) ~= P(u) + u_lo * P'(u), where P'(u) ~= K1 + 2 * K2 * u
    // (the next terms are less than 2^-8 relative to `P'(u)`, and
    // |u_lo * P'(u)| is less than 2^-59 relative to the result)
    let p_lo = p_lo + u_lo * (K1.hi() + u * (2.0 * K2.hi()));
    // exact products
    let (p_hi, p_lo) = Split::new(x).mul(p_hi, p_lo);

    // The reciprocal is the largest term
    let s1 = F64x2::add11(t_hi, p_hi);
    let s2 = F64x2::fast_add11(-r_hi, s1.hi());
    (s2.hi(), s2.lo() + (s1.lo() + ((t_lo + p_lo) - r_lo)))
}

/// Returns `(hi, lo)` such that `hi + lo ~= ln(x / c)`, where `d = x - c`
/// (exact), for positive and normal `x` and `c`, with a relative error of
/// about 2^-75 (see `ln_1p_poly_accurate` and `ln_accurate_parts`).
#[inline]
fn ln_ratio(x: f64, c: f64, d: f64) -> (f64, f64) {
    // q + q_lo = d / c
    let q = (d / c).purify();
    if q.abs() < 1.0 / 256.0 {
        // ln(x / c) = ln(1 + d / c), where the relative error of the
        // quotient is less than 2^-77
        let q_lo = div_residual(d, c, q) / c;
        ln_1p_poly_accurate(q, q_lo)
    } else {
        // w + w_lo = x / c, and the result is not close to zero
        let w = (x / c).purify();
        let w_lo = div_residual(x, c, w) / c;
        ln_accurate_sum_parts(w, w_lo)
    }
}

/// Returns `(hi, lo)` such that `hi + lo ~= y0(x)`, for `x` in the interval
/// of `poly`, which is one of the first `Y0_ANCHORED` of `Y0_TBL`.
///
/// `y0(x) = g(x) + 2/π * j0(x) * ln(x / c)`, where `poly` approximates
/// `g(x)`. The relative error is about 2^-68:
/// * `g(x)`: 2^-69 (see `ZeroPoly::eval`).
/// * `j0(x)`: 2^-69 (see `j0_parts`), in a term that is less than two
///   thirds of the result close to the zero, and less than the result
///   everywhere.
/// * Logarithm: 2^-75 (see `ln_ratio`), in the same term.
///
/// Both terms have the same sign close to the zero, where they are small.
#[inline]
fn y0_anchored(x: f64, poly: &ZeroPoly) -> (f64, f64) {
    let (g_hi, g_lo) = poly.eval(x);
    let (j_hi, j_lo) = j0_parts(x);
    let (l_hi, l_lo) = ln_ratio(x, poly.c, x - poly.c);
    let (a_hi, a_lo) = mul22(j_hi, j_lo, l_hi, l_lo);
    let (t_hi, t_lo) = mul_frac_2_pi(a_hi, a_lo);

    let s = F64x2::add11(g_hi, t_hi);
    (s.hi(), s.lo() + (g_lo + t_lo))
}

/// Returns `(hi, lo)` such that `hi + lo ~= y1(x)`, for `x` in the interval
/// of `poly`, which is one of the first `Y1_ANCHORED` of `Y1_TBL`.
///
/// `y1(x) = g(x) / x + 2/π * j1(x) * ln(x / c)`, where `poly` approximates
/// `g(x)`. The relative error is about 2^-68, like in `y0_anchored`.
#[inline]
fn y1_anchored(x: f64, poly: &ZeroPoly) -> (f64, f64) {
    let (g_hi, g_lo) = poly.eval(x);
    let (j_hi, j_lo) = j1_parts(x);
    let (l_hi, l_lo) = ln_ratio(x, poly.c, x - poly.c);
    let (a_hi, a_lo) = mul22(j_hi, j_lo, l_hi, l_lo);
    let (t_hi, t_lo) = mul_frac_2_pi(a_hi, a_lo);

    // q + q_lo = g(x) / x
    let q = (g_hi / x).purify();
    let q_lo = (div_residual(g_hi, x, q) + g_lo) / x;

    let s = F64x2::add11(q, t_hi);
    (s.hi(), s.lo() + (q_lo + t_lo))
}

/// Approximations of the modulus and the phase of the Bessel functions of an
/// order, for `x >= J0_LARGE_MIN`, with `s = 1 / x^2`:
///
/// * `M(x) ~= sqrt(2 / (π * x)) * (1 + A1 * s + ... + A7 * s^7)`
/// * `θ(x) ~= x - offset * π/128 + (B0 + B1 * s + ... + B7 * s^7) / x`
struct Large {
    /// `(2 * order + 1) * π/4`, in units of π/128
    offset: u8,
    /// `A1` and `[A2, ..., A7]`
    modulus: (F64x2, [f64; 6]),
    /// `[B0, B1]` and `[B2, ..., B7]`
    phase: ([F64x2; 2], [f64; 6]),
    /// More accurate approximation of the phase, with more terms and fixed
    /// point coefficients (see `phase_precise`)
    phase_fixed: &'static [i128; 17],
}

static LARGE0: Large = Large {
    offset: 32,
    modulus: {
        // GENERATE: bessel_modulus_poly 0 F64x2:1,f64 7 47.875
        const K1: F64x2 = F64x2::from_bits(0xBFB0000000000000, 0x3BE9CCFD5C903A9D); // -6.249999999999999995629184463437e-2
        const K2: f64 = f64::from_bits(0x3FBA7FFFFFFFFF47); // 1.0351562499999743e-1
        const K3: f64 = f64::from_bits(0xBFE15EFFFFF8FDEC); // -5.428466796365092e-1
        const K4: f64 = f64::from_bits(0x401765115FBAA19F); // 5.848699088827202e0
        const K5: f64 = f64::from_bits(0xC05AB8996905515F); // -1.068843634178097e2
        const K6: f64 = f64::from_bits(0x40A722B8E73562AF); // 2.961361138981148e3
        const K7: f64 = f64::from_bits(0xC0FA07AB3697EE39); // -1.0661870082848602e5

        (K1, [K2, K3, K4, K5, K6, K7])
    },
    phase: {
        // GENERATE: bessel_phase_poly 0 F64x2:2,f64 8 47.875
        const K0: F64x2 = F64x2::from_bits(0xBFC0000000000000, 0x0000000000000000); // -1.250000000000000000000000000000e-1
        const K1: F64x2 = F64x2::from_bits(0x3FB0AAAAAAAAAAAB, 0xBC58AA5E839295AA); // 6.510416666666666594411718652823e-2
        const K2: f64 = f64::from_bits(0xBFCAD33333332E72); // -2.0957031249996622e-1
        const K3: f64 = f64::from_bits(0x3FFA3584922185BE); // 1.6380658825149648e0
        const K4: f64 = f64::from_bits(0xC03779A1A7340FD8); // -2.3475122881136343e1
        const K5: f64 = f64::from_bits(0x4080BCF1705F2CA8); // 5.356178901134017e2
        const K6: f64 = f64::from_bits(0xC0D15C88E1EF5EC7); // -1.7778138789980465e4
        const K7: f64 = f64::from_bits(0x41266A8698D89156); // 7.345312985272806e5

        ([K0, K1], [K2, K3, K4, K5, K6, K7])
    },
    phase_fixed: &PHASE0_FIXED,
};

static LARGE1: Large = Large {
    offset: 96,
    modulus: {
        // GENERATE: bessel_modulus_poly 1 F64x2:1,f64 7 47.875
        const K1: F64x2 = F64x2::from_bits(0x3FC8000000000000, 0xBBED712BCCB5407F); // 1.874999999999999999501231855913e-1
        const K2: f64 = f64::from_bits(0xBFC8BFFFFFFFFF96); // -1.9335937499999706e-1
        const K3: f64 = f64::from_bits(0x3FE9C4FFFFF800F6); // 8.052978515043197e-1
        const K4: f64 = f64::from_bits(0xC01EF5B65B2E9174); // -7.739953446118637e0
        const K5: f64 = f64::from_bits(0x4060984A26727C96); // 1.3275905153618504e2
        const K6: f64 = f64::from_bits(0xC0AB9F231AF9EBF9); // -3.535568565187507e3
        const K7: f64 = f64::from_bits(0x40FE4B40D7AD17F1); // 1.2408405265530922e5

        (K1, [K2, K3, K4, K5, K6, K7])
    },
    phase: {
        // GENERATE: bessel_phase_poly 1 F64x2:2,f64 8 47.875
        const K0: F64x2 = F64x2::from_bits(0x3FD8000000000000, 0x0000000000000000); // 3.750000000000000000000000000000e-1
        const K1: F64x2 = F64x2::from_bits(0xBFC5000000000000, 0x3C2E50236C08B8CC); // -1.640624999999999991783633909642e-1
        const K2: f64 = f64::from_bits(0x3FD7BCCCCCCCCA19); // 3.708984374999616e-1
        const K3: f64 = f64::from_bits(0xC002F486DB5730B2); // -2.369397843931517e0
        const K4: f64 = f64::from_bits(0x403E9FBEE321DE15); // 3.0624006458061633e1
        const K5: f64 = f64::from_bits(0xC0849946A6CB1284); // -6.591594978203798e2
        const K6: f64 = f64::from_bits(0x40D49846564D78CD); // 2.1089099017494173e4
        const K7: f64 = f64::from_bits(0xC129F81316C7C118); // -8.509535444927542e5

        ([K0, K1], [K2, K3, K4, K5, K6, K7])
    },
    phase_fixed: &PHASE1_FIXED,
};

/// Number of fractional bits of the coefficients of `Large::phase_fixed`.
const PHASE_FIXED_BITS: u32 = 123;

// GENERATE: bessel_phase_fixed 0 17 47.875 123
// PHASE0_FIXED[k] = 128/π * K[k] * 2^(-10 * k), in units of 2^-123, where
// θ(x) ~= x - π/4 + (K[0] + K[1] / x^2 + ... + K[16] / x^32) / x, for x >= 47.875
static PHASE0_FIXED: [i128; 17] = [
    -0x28BE60DB9391054A7F09D5F47D4D3770, // -5.092958178940651e0
    0x00054E1F4741E18CB05DDFF28680504B,  // 2.5904124852196508e-3
    -0x00000444F2E63556789454C301716973, // -8.143103310470757e-6
    0x0000000857B4B39935B11A18215F53E0,  // 6.215722141366143e-8
    -0x000000001DE3B3B6CCD130A521CAF511, // -8.698978044842151e-10
    0x0000000000AA7FEA8E825AB69B37921E,  // 1.938356863329894e-11
    -0x0000000000058B7213F30E9BB9FD2D8C, // -6.303604824361751e-13
    0x0000000000003F787BAC66183BDF85D2,  // 2.8186625081449692e-14
    -0x00000000000003BC7F85B3625DFE5DBB, // -1.6592597675368214e-15
    0x0000000000000047B68DED4E12330E22,  // 1.2440240354123772e-16
    -0x0000000000000006AC071E1A433B11C6, // -1.1574046592274656e-17
    0x0000000000000000C1186652BB35A8A8,  // 1.3084647252387914e-18
    -0x00000000000000001A0FD13E5C99DD1C, // -1.766015350041221e-19
    0x0000000000000000041F0F0B176D8B00,  // 2.792717292326608e-20
    -0x000000000000000000BB7102DEF7EACB, // -4.961533871245345e-21
    0x00000000000000000020BB8249F8F896,  // 8.664209173626735e-22
    -0x00000000000000000003DE4FCB0D54F6, // -1.0239582900902308e-22
];

// GENERATE: bessel_phase_fixed 1 17 47.875 123
// PHASE1_FIXED[k] = 128/π * K[k] * 2^(-10 * k), in units of 2^-123, where
// θ(x) ~= x - 3 * π/4 + (K[0] + K[1] / x^2 + ... + K[16] / x^32) / x, for x >= 47.875
static PHASE1_FIXED: [i128; 17] = [
    0x7A3B2292BAB30FDF7D1D81DD77E7A651,  // 1.5278874536821952e1
    -0x000D5E77C80C6B95BC71AF3A34391CE2, // -6.527839462753519e-3
    0x0000078E4E3651FDB80E63E8BFC05CAE,  // 1.4411699148726904e-5
    -0x0000000C1136E952BB34C06AE51CD085, // -8.990797498638916e-8
    0x0000000026FDE040EC3D7D78D6A1834D,  // 1.1348079159055844e-9
    -0x0000000000D1D33E8887C3EFD6412435, // -2.385436038511715e-11
    0x0000000000069390E3FF8458820E7679,  // 7.476534853335063e-13
    -0x00000000000049632DC9F1EE03BAB91B, // -3.2590560220477936e-14
    0x000000000000043D7149DFF37B7A26CF,  // 1.882942644669066e-15
    -0x00000000000000503845035FF7BADFC8, // -1.391591756022608e-16
    0x000000000000000760FA9E6BB0C818B7,  // 1.2800219460542342e-17
    -0x0000000000000000D3949CC65BA76937, // -1.4337253524844224e-18
    0x00000000000000001C56812E0120DD59,  // 1.9202513809911892e-19
    -0x00000000000000000473DC1CE8B14AFB, // -3.017183810941201e-20
    0x000000000000000000C975EE69F2ABCD,  // 5.332619510990495e-21
    -0x000000000000000000230CA47CBC395F, // -9.277494928482351e-22
    0x0000000000000000000421F7A4DC47D4,  // 1.0939126104304747e-22
];

/// Returns `(hi, lo)` such that `hi + lo ~= M(x) * sin(θ(x) + quarter * π/128)`,
/// for `x >= J0_LARGE_MIN`, where `M` and `θ` are the modulus and the phase
/// of `large` (see `Large::eval`): `quarter` is 0 for the function of the
/// second kind and 64 for the one of the first kind.
///
/// The relative error is about 2^-65:
/// * Sine: 2^-65, from the terms that `Reduced::sin_k` evaluates with plain
///   `f64` arithmetic (the largest one is less than 2^-13.7 relative to the
///   result and has three rounding errors).
/// * Phase: less than 2^-71 relative to its distance to the zero.
/// * Modulus: 2^-74.
#[inline(never)]
fn large_eval(x: f64, large: &Large, quarter: u8) -> (f64, f64) {
    large_eval_arg(&LargeArg::new(x), large, quarter)
}

/// Returns `(large_eval(x, &LARGE0, quarter), large_eval(x, &LARGE1, quarter))`
/// (see `LargeArg`).
#[inline(never)]
fn large_eval_both(x: f64, quarter: u8) -> ((f64, f64), (f64, f64)) {
    let arg = LargeArg::new(x);
    (
        large_eval_arg(&arg, &LARGE0, quarter),
        large_eval_arg(&arg, &LARGE1, quarter),
    )
}

/// Like `large_eval`, from the parts that do not depend on the order.
#[inline(always)]
fn large_eval_arg(arg: &LargeArg, large: &Large, quarter: u8) -> (f64, f64) {
    let (m_hi, m_lo, mut n, mut b) = large.eval(arg);
    if n.wrapping_add(quarter) & 127 == 0 && b.bh.abs() < NEAR_ZERO {
        (n, b) = phase_precise(arg.x, large);
    }
    // The low part of the sine is not small enough for the product
    let t = b.sin_k(n.wrapping_add(quarter)).normalize();
    mul22(m_hi, m_lo, t.hi(), t.lo())
}

/// Returns `j0(x)`, for `J0_LARGE_MIN <= x <= 2^128`, with the relative error
/// of `Large::eval_short`.
#[inline]
pub(crate) fn j0_large_short(x: f64) -> f64 {
    LARGE0.eval_short(x, 64)
}

/// Returns `j1(x)`, for `J0_LARGE_MIN <= x <= 2^128`, with the relative error
/// of `Large::eval_short`.
#[inline]
pub(crate) fn j1_large_short(x: f64) -> f64 {
    LARGE1.eval_short(x, 64)
}

/// Returns `y0(x)`, for `J0_LARGE_MIN <= x <= 2^128`, with the relative error
/// of `Large::eval_short`.
#[inline]
pub(crate) fn y0_large_short(x: f64) -> f64 {
    LARGE0.eval_short(x, 0)
}

/// Returns `y1(x)`, for `J0_LARGE_MIN <= x <= 2^128`, with the relative error
/// of `Large::eval_short`.
#[inline]
pub(crate) fn y1_large_short(x: f64) -> f64 {
    LARGE1.eval_short(x, 0)
}

/// Returns `(y0(x), y1(x))`, for `J0_LARGE_MIN <= x <= 2^128`, like
/// `y0_large_short` and `y1_large_short`, sharing the parts that do not
/// depend on the order (see `ShortArg`).
#[inline]
pub(crate) fn y0_y1_large_short(x: f64) -> (f64, f64) {
    let arg = ShortArg::new(x);
    (
        LARGE0.eval_short_arg(&arg, 0),
        LARGE1.eval_short_arg(&arg, 0),
    )
}

/// The parts of `Large::eval_short` that do not depend on the order, so the
/// functions of orders 0 and 1 can share them (like `LargeArg`).
#[derive(Clone, Copy)]
struct ShortArg {
    x_exp: i16,
    /// `x = n * π/128 + b`
    n: u8,
    b: Reduced,
    /// `1 / x = t_hi + t_lo`
    t_hi: f64,
    t_lo: f64,
    /// `sqrt(2 / (π * x))`
    r: f64,
    /// `1 / x^2`
    s: f64,
}

impl ShortArg {
    #[inline(always)]
    fn new(x: f64) -> Self {
        // x = n * π/128 + b
        let (n, b) = reduce_rad(x);

        // 1 / x = t_hi + t_lo
        let (t_hi, t_lo) = recip_parts(x);
        // sqrt(2 / (π * x))
        let (r, _) = fast_sqrt(FRAC_2_PI.hi() * t_hi);

        Self {
            x_exp: x.exponent(),
            n,
            b,
            t_hi,
            t_lo,
            r,
            s: t_hi * t_hi,
        }
    }
}

/// The parts of `Large::eval` that do not depend on the order, for
/// `x >= J0_LARGE_MIN`, so the functions of orders 0 and 1 can share them
/// (see `jn_forward` and `yn_finite`).
struct LargeArg {
    x: f64,
    x_exp: i16,
    /// `x = n * π/128 + b`
    n: u8,
    b: Reduced,
    /// `1 / x = (t_hi + t_lo) * m_scale^2`
    t_hi: f64,
    t_lo: f64,
    /// `sqrt(2 / (π * x)) = (r_hi + r_lo) * m_scale`
    r_hi: f64,
    r_lo: f64,
    m_scale: f64,
}

impl LargeArg {
    #[inline(always)]
    fn new(x: f64) -> Self {
        // x = n * π/128 + b
        let (n, b) = reduce_rad(x);
        let x_exp = x.exponent();

        // 1 / x = (t_hi + t_lo) * m_scale^2. When `x` is huge, it is scaled
        // so its reciprocal and its products are not subnormal.
        let (xs, m_scale) = if x_exp < LARGE_SCALED_MIN_EXP {
            (x, 1.0)
        } else {
            (x * f64::exp2i_fast(-512), f64::exp2i_fast(-256))
        };
        let (t_hi, t_lo) = recip_parts(xs);

        // sqrt(2 / (π * x)) = (r_hi + r_lo) * m_scale
        let (w_hi, w_lo) = mul_frac_2_pi(t_hi, t_lo);
        let w = F64x2::fast_add11(w_hi, w_lo);
        let (r_hi, r_lo) = sqrt_parts(w.hi(), w.lo());

        Self {
            x,
            x_exp,
            n,
            b,
            t_hi,
            t_lo,
            r_hi,
            r_lo,
            m_scale,
        }
    }
}

impl Large {
    /// Returns `M(x) * sin(θ(x) + quarter * π/128)`, like `large_eval`, for
    /// `J0_LARGE_MIN <= x <= 2^128`.
    ///
    /// It is meant for the `f32` functions: the polynomials and the modulus
    /// are evaluated with plain `f64` arithmetic, so the relative error of
    /// the modulus is about 2^-50, but the phase is still calculated as a
    /// sum of two `f64`. Its absolute error is the one of its polynomial
    /// (less than 2^-82.5) plus the one of the evaluation of the terms after
    /// the first one (less than 2^-11 relative to it), which is about 2^-64
    /// relative to the last term of the phase (less than `3 / (8 * x)`).
    ///
    /// That error is less than 2^-46 relative to the distance from the phase
    /// to a zero of the function with every `f32` value (found by exhaustive
    /// search: the closest ones to a zero where that term is not much
    /// smaller are at 2^-25.6), so the relative error of the result is about
    /// 2^-46 there, and 2^-50 away from the zeros.
    #[inline]
    fn eval_short(&self, x: f64, quarter: u8) -> f64 {
        self.eval_short_arg(&ShortArg::new(x), quarter)
    }

    /// Like `eval_short`, from the parts that do not depend on the order.
    #[inline(always)]
    fn eval_short_arg(&self, arg: &ShortArg, quarter: u8) -> f64 {
        let ShortArg {
            x_exp,
            n,
            b,
            t_hi,
            t_lo,
            r,
            s,
        } = *arg;
        let n = n.wrapping_sub(self.offset);

        let ([b0, b1], [b2, b3, b4, b5, b6, b7]) = self.phase;
        let (a, q_lo) = if x_exp < LARGE_POLY_MAX_EXP {
            let s2 = s * s;
            let s4 = s2 * s2;
            let (a1, [a2, a3, a4, a5, a6, a7]) = self.modulus;
            let wa = ((a2 + s * a3) + s2 * (a4 + s * a5)) + s4 * (a6 + s * a7);
            let wb = ((b2 + s * b3) + s2 * (b4 + s * b5)) + s4 * (b6 + s * b7);
            (1.0 + s * (a1.hi() + s * wa), s * (b1.hi() + s * wb))
        } else {
            // The terms after the first one of each polynomial are
            // negligible, and skipping them avoids the powers of `1 / x`
            // that are subnormal
            (1.0, 0.0)
        };

        // The last term of the phase is `(B0 + q_lo) / x`, where `B0` has
        // few bits, so its products are exact
        let (p_hi, p_lo) = mul22(t_hi, t_lo, b0.hi(), q_lo);
        let (n, b) = add_phase(n, b, p_hi, p_lo);
        (r * a) * b.sin_k(n.wrapping_add(quarter)).to_f64()
    }

    /// Returns `(m_hi, m_lo, n, b)` such that `m_hi + m_lo ~= M(x)` and
    /// `θ(x) ~= n * π/128 + b` (modulo `2 * π`), for `x >= J0_LARGE_MIN`,
    /// from the parts that do not depend on the order (see `LargeArg`).
    ///
    /// The relative error of the modulus is about 2^-74:
    /// * Polynomial approximation: less than 2^-78 (2^-81 without the
    ///   rounding of the coefficients).
    /// * Reciprocal: 2^-78 (half the error of `recip_parts`).
    /// * Evaluation of the polynomial and the products: ~2^-75.
    ///
    /// The absolute error of the phase is about 2^-81:
    /// * Polynomial approximation: less than 2^-82.5 (2^-85 without the
    ///   rounding of the coefficients).
    /// * Reciprocal: 2^-77 relative to the last term, which is less than
    ///   2^-7.
    /// * Evaluation of the terms of degree 2 and higher of the polynomial
    ///   (less than 2^-22 relative to that term): ~2^-81.
    /// * Reduction of `x`: less than 2^-100.
    #[inline(always)]
    fn eval(&self, arg: &LargeArg) -> (f64, f64, u8, Reduced) {
        let LargeArg {
            x_exp,
            n,
            b,
            t_hi,
            t_lo,
            r_hi,
            r_lo,
            m_scale,
            ..
        } = *arg;
        let n = n.wrapping_sub(self.offset);

        if x_exp >= LARGE_POLY_MAX_EXP {
            let m = (r_hi * m_scale, r_lo * m_scale);
            if x_exp >= LARGE_PHASE_MAX_EXP {
                return (m.0, m.1, n, b);
            }
            // θ(x) ~= x - offset * π/128 + B0 / x, where `1 / x` is not scaled
            // (see `LARGE_SCALED_MIN_EXP`)
            let ([b0, _], _) = self.phase;
            let (p_hi, p_lo) = mul22(t_hi, t_lo, b0.hi(), b0.lo());
            let (n, b) = add_phase(n, b, p_hi, p_lo);
            return (m.0, m.1, n, b);
        }

        // The parts of `arg` are not scaled here (see `LARGE_SCALED_MIN_EXP`)

        // s + s_lo = 1 / x^2, where `s` is `t_hi^2` rounded
        let (s, se) = square_parts(t_hi);
        let s_lo = se + 2.0 * t_hi * t_lo;
        let s2 = s * s;
        let s4 = s2 * s2;
        let ss = Split::new(s);

        // The terms of degree 2 and higher of both polynomials are
        // evaluated with plain `f64` arithmetic, and the next steps use
        // exact products
        let (a1, [a2, a3, a4, a5, a6, a7]) = self.modulus;
        let wa = ((a2 + s * a3) + s2 * (a4 + s * a5)) + s4 * (a6 + s * a7);
        let (h, l) = ss.mul_add(wa, 0.0, a1);
        let (a_hi, a_lo) = ss.mul_add(h, l, F64x2::new1(1.0));
        // A(s + s_lo) ~= A(s) + s_lo * A1 (the next terms are negligible)
        let a_lo = a_lo + s_lo * a1.hi();
        let (m_hi, m_lo) = mul22(r_hi, r_lo, a_hi, a_lo);

        let ([b0, b1], [b2, b3, b4, b5, b6, b7]) = self.phase;
        let wb = ((b2 + s * b3) + s2 * (b4 + s * b5)) + s4 * (b6 + s * b7);
        let (h, l) = ss.mul_add(wb, 0.0, b1);
        let (q_hi, q_lo) = ss.mul_add(h, l, b0);
        // B(s + s_lo) ~= B(s) + s_lo * B1 (the next terms are negligible)
        let q_lo = q_lo + s_lo * b1.hi();
        let (p_hi, p_lo) = mul22(t_hi, t_lo, q_hi, q_lo);
        let (n, b) = add_phase(n, b, p_hi, p_lo);

        (m_hi, m_lo, n, b)
    }
}

/// Returns `(k, c)` such that `k * π/128 + c = n * π/128 + b + p_hi + p_lo`,
/// with `|c| <= ~π/256`, when `|b + p_hi + p_lo|` is at most about π/128 (as
/// when `|b|` and `|p_hi + p_lo|` are at most about π/256): one addition or
/// subtraction of π/128 is enough then, and it is exact (by Sterbenz lemma).
#[inline]
fn add_phase(n: u8, b: Reduced, p_hi: f64, p_lo: f64) -> (u8, Reduced) {
    let s = F64x2::add11(b.bh, p_hi);
    let s = F64x2::add11(s.hi(), s.lo() + (b.bl + p_lo));
    let (hi, lo) = (s.hi(), s.lo());
    if hi.abs() <= 0.5 * FRAC_PI_128_R0 {
        (n, Reduced { bh: hi, bl: lo })
    } else if hi.is_sign_negative() {
        // The sum is exact (by Sterbenz lemma)
        let bh = hi + FRAC_PI_128_R0;
        let bl = lo + FRAC_PI_128_R1;
        (n.wrapping_sub(1), Reduced { bh, bl })
    } else {
        // The subtraction is exact (by Sterbenz lemma)
        let bh = hi - FRAC_PI_128_R0;
        let bl = lo - FRAC_PI_128_R1;
        (n.wrapping_add(1), Reduced { bh, bl })
    }
}

/// Returns `(n, b)` such that `θ(x) ~= n * π/128 + b` (modulo `2 * π`), with
/// `|b| <= π/256`, for `x >= J0_LARGE_MIN`, where `θ` is the phase of
/// `large`.
///
/// The phase is calculated in units of π/128 with integer arithmetic:
///
/// `θ(x) * 128/π = x * 128/π - offset + y * 2^-e * sum(D[k] * σ^k)`
///
/// where `2^e <= x < 2^(e + 1)`, `y = 2^e / x`, `σ = 2^10 / x^2` and `D` is
/// `large.phase_fixed`. The absolute error of `b` is about 2^-130 (when it
/// is small):
/// * Polynomial approximation: 2^-139.
/// * `x * 128/π`: 2^-128 (see `mul_128_pi`), which is 2^-133.3 radians.
/// * Evaluation of the polynomial: each step has an error of up to 4 units
///   of 2^-123, and the next ones multiply it by `σ < 0.45`, so they add less
///   than 8 units, which are 2^-125 after the product by `y * 2^-e` (less than
///   2^-5.5), or 2^-130.3 radians.
#[cold]
fn phase_precise(x: f64, large: &Large) -> (u8, Reduced) {
    /// Returns the high 128 bits of `a * b`, with an error of up to 3 (it
    /// ignores the product of the low halves and the carries of the others).
    #[inline]
    fn mul_hi(a: u128, b: u128) -> u128 {
        let (a_hi, a_lo) = (a >> 64, a & u128::from(u64::MAX));
        let (b_hi, b_lo) = (b >> 64, b & u128::from(u64::MAX));
        a_hi * b_hi + ((a_hi * b_lo) >> 64) + ((a_lo * b_hi) >> 64)
    }

    /// Like `mul_hi`, with a signed `a`.
    #[inline]
    fn mul_hi_signed(a: i128, b: u128) -> i128 {
        let p = mul_hi(a.unsigned_abs(), b) as i128;
        if a < 0 { -p } else { p }
    }

    // x * 128/π = n + frac * 2^-128 (modulo 256)
    let (n, frac) = mul_128_pi(x);

    // x = m * 2^(e - 52), with e >= 5
    let e = x.exponent() as u32;
    let m = u128::from(x.mant());
    // yq = y * 2^127 = 2^179 / m, where 1/2 < y <= 1
    let one = 1u128 << 115;
    let yq = ((one / m) << 64) | (((one % m) << 64) / m);
    // sq = σ * 2^128 = y^2 * 2^(138 - 2 * e), where σ < 0.45
    let y2 = mul_hi(yq, yq);
    let sq = if e == 5 {
        y2 << 2
    } else {
        y2.checked_shr(2 * e - 12).unwrap_or(0)
    };

    // acc = sum(D[k] * σ^k), in units of 2^-PHASE_FIXED_BITS
    let (last, others) = large.phase_fixed.split_last().unwrap();
    let acc = others
        .iter()
        .rev()
        .fold(*last, |acc, &d| d + mul_hi_signed(acc, sq));

    // p = y * 2^-e * acc, in units of 2^-128: |p| < 2^-1.6 and the product
    // of `acc` and `yq` is in units of 2^(-PHASE_FIXED_BITS + 1)
    let p = mul_hi_signed(acc, yq);
    let shift = (128 - PHASE_FIXED_BITS + 1) as i32 - e as i32;
    let p = if shift >= 0 {
        p << shift
    } else {
        p.checked_shr(shift.unsigned_abs()).unwrap_or(0)
    };

    // n + frac * 2^-128 + p * 2^-128
    let (frac, carry) = frac.overflowing_add(p as u128);
    let n = if p >= 0 {
        n.wrapping_add(u8::from(carry))
    } else {
        n.wrapping_sub(u8::from(!carry))
    };
    reduced_from_fixed(n.wrapping_sub(large.offset), frac)
}

/// `jn(n, x)` and `yn(n, x)` are calculated from their modulus and phase for
/// `x >= n^2 * ORDER_LARGE_FACTOR` (see `order_large_eval`).
const ORDER_LARGE_FACTOR: f64 = 65536.0;

/// With large arguments and other orders, the terms after the first one of
/// the modulus and the phase are only evaluated when the exponent of `x` is
/// less than this: the largest one, `(4 * n^2 - 1) / (8 * x) < 2^61 / x`, is
/// less than 2^-139 otherwise. Below it, the powers of `1 / x` are not
/// subnormal.
const ORDER_LARGE_POLY_MAX_EXP: i16 = 200;

// `order_large_eval` uses `1 / x` from `LargeArg` below
// `ORDER_LARGE_POLY_MAX_EXP`, where it is not scaled
const _: () = assert!(ORDER_LARGE_POLY_MAX_EXP <= LARGE_SCALED_MIN_EXP);

/// Returns `(hi, lo)` such that `hi + lo ~= M(x) * sin(θ(x) + quarter * π/128)`,
/// like `large_eval`, where `M` and `θ` are the modulus and the phase of the
/// functions of order `n` (see `order_large_parts`), for `n >= 2` and
/// `x >= n^2 * ORDER_LARGE_FACTOR`.
///
/// The relative error is about 2^-65, the one of the sine (see
/// `large_eval`). Unlike with orders 0 and 1, the phase is not calculated
/// again when it is close to a zero of the function (see `phase_precise`),
/// so the error of the results close to a zero is about 2^-88 times the
/// modulus (the error of the phase).
#[inline(never)]
fn order_large_eval(n: u32, x: f64, quarter: u8) -> (f64, f64) {
    let (m_hi, m_lo, k, b) = order_large_parts(n, x);
    // The low part of the sine is not small enough for the product
    let t = b.sin_k(k.wrapping_add(quarter)).normalize();
    mul22(m_hi, m_lo, t.hi(), t.lo())
}

/// Returns `(m_hi, m_lo, k, b)` such that `m_hi + m_lo ~= M(x)` and
/// `θ(x) ~= k * π/128 + b` (modulo `2 * π`), where `M` and `θ` are the
/// modulus and the phase of the functions of order `n`, for `n >= 2` and
/// `x >= n^2 * ORDER_LARGE_FACTOR`.
///
/// They are the first terms of their asymptotic series, with `μ = 4 * n^2`
/// and `s = 1 / x^2`:
///
/// * `M(x) ~= sqrt(2 / (π * x)) * (1 + A1 * s + A2 * s^2)`, where
///   `A1 = (μ - 1) / 16` and `A2 = (μ - 1) * (5 * μ - 53) / 512`
/// * `θ(x) ~= x - (2 * n + 1) * π/4 + (B0 + B1 * s) / x`, where
///   `B0 = (μ - 1) / 8` and `B1 = (μ - 1) * (μ - 25) / 384`
///
/// which are the first coefficients of the polynomials of `Large` with
/// `μ = 0` and `μ = 4` (orders 0 and 1). With these arguments, the next
/// terms are less than about `n^6 / (8 * x^6)` (modulus, relative to it) and
/// `n^6 / (80 * x^5)` (phase), so the truncation errors are less than 2^-106
/// and 2^-89 (with `n = 2`, and much less with larger orders).
///
/// The relative error of the modulus is about 2^-76, the one of
/// `sqrt(2 / (π * x))` in `LargeArg` (half the errors of `1 / x` and its
/// product by 2/π). The absolute error of the phase is about 2^-88:
/// * Truncation: less than 2^-89.
/// * `B0 / x`, which is less than 2^-17: 2^-75 relative to it (the error of
///   `mul22`, plus the one of `recip_parts`).
/// * `B1 * s / x`, which is less than 2^-55: 2^-50 relative to it.
/// * Reduction of `x`: less than 2^-100.
#[inline(always)]
fn order_large_parts(n: u32, x: f64) -> (f64, f64, u8, Reduced) {
    let LargeArg {
        x_exp,
        n: k,
        b,
        t_hi,
        t_lo,
        r_hi,
        r_lo,
        m_scale,
        ..
    } = LargeArg::new(x);
    // (2 * n + 1) * π/4 = (64 * n + 32) * π/128 (modulo 2 * π)
    let k = k.wrapping_sub((n as u8).wrapping_mul(64).wrapping_add(32));

    if x_exp >= ORDER_LARGE_POLY_MAX_EXP {
        (r_hi * m_scale, r_lo * m_scale, k, b)
    } else {
        // `1 / x` and `sqrt(2 / (π * x))` are not scaled here (see
        // `ORDER_LARGE_POLY_MAX_EXP`)

        // n^2 = nn.hi() + nn.lo() exactly, where `nn.lo()` is zero or an
        // integer, so B0 = n^2 / 2 - 1/8 = c.hi() + c.lo() exactly
        let nf = f64::from(n);
        let nn = F64x2::mul11(nf, nf);
        let c = F64x2::add11(0.5 * nn.hi(), 0.5 * nn.lo() - 0.125);
        let (p_hi, p_lo) = mul22(t_hi, t_lo, c.hi(), c.lo());

        // The other terms are evaluated with plain `f64` arithmetic: they
        // are less than 2^-55 (phase) and 2^-36 (modulus, relative to it)
        let mu = 4.0 * nn.hi();
        let s = t_hi * t_hi;
        let q = (((mu - 1.0) * (mu - 25.0)) * (1.0 / 384.0)) * s * t_hi;
        let (k, b) = add_phase(k, b, p_hi, p_lo + q);
        let a1 = (mu - 1.0) * (1.0 / 16.0);
        let a2 = ((mu - 1.0) * (5.0 * mu - 53.0)) * (1.0 / 512.0);
        let a = s * (a1 + s * a2);
        (r_hi, r_lo + r_hi * a, k, b)
    }
}

/// `jn(n, x)` and `yn(n, x)` are calculated from Debye's asymptotic
/// expansions for `n >= ORDER_DEBYE_MIN` and `x >= DEBYE_BELOW_MIN * n`
/// (unless `x` is large enough for `order_large_eval`), which take a time that
/// does not depend on the order (see `jn_debye` and `yn_debye`, which use
/// Olver's expansions close to `x = n`).
///
/// From this order, their truncation errors are small enough (see
/// `DEBYE_DIST_MIN` and `olver`), and they are much faster than the
/// recurrences.
pub(crate) const ORDER_DEBYE_MIN: u32 = 1024;

/// `debye_below` is only used for `x >= DEBYE_BELOW_MIN * n`, so
/// `w = sqrt(1 - (x / n)^2) < 0.933`, where `atanh_precise` is accurate.
///
/// With `n >= ORDER_DEBYE_MIN`, the results with smaller arguments are zero
/// (`jn`) or infinite (`yn`), which `jn_underflows` and `yn_overflows` detect
/// (see `test_debye_below_min`), so no order from it needs the recurrences.
/// With smaller orders, they are calculated with the recurrences (see
/// `jn_backward` and `yn_finite`).
pub(crate) const DEBYE_BELOW_MIN: f64 = 0.36;

// `atanh_precise` is accurate up to the last point of `ATANH_TBL`
const _: () = {
    let w_max = (ATANH_TBL.len() - 1) as f64 / 128.0;
    assert!(1.0 - DEBYE_BELOW_MIN * DEBYE_BELOW_MIN < w_max * w_max);
};

/// The Debye expansions are only used when `R^3 / n^2 >= DEBYE_DIST_MIN`,
/// where `R = sqrt(|n^2 - x^2|)`, so their variable `a = n^2 / R^3` is at
/// most `1 / DEBYE_DIST_MIN` (see `debye_poly`). The other one, `b = 1 / R`,
/// is less than `a` when `R < n`, and less than `1 / n` otherwise.
///
/// With 14 terms, the first one that is left out is less than 2^-73 then
/// with `n >= ORDER_DEBYE_MIN` (2^-75.6 with `n >= 4096`, and less with
/// smaller `a`), which is about the size of the truncation errors with
/// `R^3 / n^2 = DEBYE_DIST_MIN`.
const DEBYE_DIST_MIN: f64 = 256.0;

/// With the orders from `ORDER_DEBYE_FAR_MIN` below `ORDER_DEBYE_MIN`,
/// `jn(n, x)` and `yn(n, x)` are calculated by `debye_above` when
/// `R = sqrt(x^2 - n^2)` is at least this (see `debye_far`), instead of the
/// recurrence of `jn_forward`, which takes a time proportional to the order.
///
/// Then the variables of `debye_poly`, `a = n^2 / R^3` and `b = 1 / R`, are
/// less than with `n = ORDER_DEBYE_MIN` at the limit of `debye_valid` (see
/// below), and the terms of degree `k` of the expansion are at most
/// `sum(|C[k][i]| * a^(k - i) * b^i)`. So the truncation error is as small
/// (the first term that is left out is less than 2^-73.2), and so are the
/// others (see `debye_above_parts`).
const DEBYE_FAR_R_MIN: f64 = 646.0;

// `R >= DEBYE_FAR_R_MIN` implies that `a` and `b` are less than with
// `R^3 = DEBYE_DIST_MIN * ORDER_DEBYE_MIN^2` and `n = ORDER_DEBYE_MIN`
const _: () = {
    let n = ORDER_DEBYE_MIN as f64;
    let r = DEBYE_FAR_R_MIN;
    assert!((r * r) * r >= DEBYE_DIST_MIN * (n * n));
};

/// Smallest order with which `debye_above` is used far beyond the order (see
/// `DEBYE_FAR_R_MIN`). With smaller ones, the recurrence of `jn_forward` is
/// about as fast or faster.
const ORDER_DEBYE_FAR_MIN: u32 = 20;

/// Number of terms of the Debye expansions (see `DEBYE_POLY`).
const DEBYE_TERMS: usize = 14;

// GENERATE: consts F64x2 LN_2PI LN_FRAC_PI_2 FRAC_1_3 FRAC_1_5
const LN_2PI: F64x2 = F64x2::from_bits(0x3FFD67F1C864BEB5, 0xBC965B5A1B7FF5DF); // 1.837877066409345483560659472811e0
const LN_FRAC_PI_2: F64x2 = F64x2::from_bits(0x3FDCE6BB25AA1316, 0xBC6DCD49C8E5AFF6); // 4.515827052894548647261952298949e-1
const FRAC_1_3: F64x2 = F64x2::from_bits(0x3FD5555555555555, 0x3C75555555555555); // 3.333333333333333333333333333333e-1
const FRAC_1_5: F64x2 = F64x2::from_bits(0x3FC999999999999A, 0xBC6999999999999A); // 2.000000000000000000000000000000e-1

// GENERATE: atan_table F64x2 7
// ATAN_TBL[i] = atan(i / 128)
static ATAN_TBL: [F64x2; 129] = [
    F64x2::from_bits(0x0000000000000000, 0x0000000000000000), // 0
    F64x2::from_bits(0x3F7FFFD555BBBA97, 0x3C068062351FBBE6), // 7.812341060101111296463391842199e-3
    F64x2::from_bits(0x3F8FFF555BBB729B, 0xBC2220C39D4DFF50), // 1.562372862047683080280152125657e-2
    F64x2::from_bits(0x3F97FEE0184A5C36, 0xBC343189FC0A354B), // 2.343320987946758445107229863742e-2
    F64x2::from_bits(0x3F9FFD55BBA97625, 0xBC35EC431444912C), // 3.123983343026827625371174489249e-2
    F64x2::from_bits(0x3FA3FD65F169C9D9, 0x3C27230A716461B5), // 3.904264995516699345275445132693e-2
    F64x2::from_bits(0x3FA7FB818430DA2A, 0xBC086EF8F794F105), // 4.684071291596965375222376000186e-2
    F64x2::from_bits(0x3FABF8DDF139C444, 0xBC489FE34B2A7FA8), // 5.463307923935947507138879374733e-2
    F64x2::from_bits(0x3FAFF55BB72CFDEA, 0xBC3C934D86D23F1D), // 6.241880999595734847397911298551e-2
    F64x2::from_bits(0x3FB1F86DBF082D59, 0xBC4095DC7732EF81), // 7.019697107187051852187630470350e-2
    F64x2::from_bits(0x3FB3F59F0E7C559D, 0x3C5AC4CE285DF847), // 7.796663383154230656332864878103e-2
    F64x2::from_bits(0x3FB5F2324FD2D7B2, 0x3C58A8DA4401318E), // 8.572687577074481459363165573197e-2
    F64x2::from_bits(0x3FB7EE182602F10F, 0xBC5CFB654C0C3D98), // 9.347678115858946350452719331206e-2
    F64x2::from_bits(0x3FB9E94153CFDCF1, 0x3C5A332E1D69C47E), // 1.012154416674666741662006965791e-1
    F64x2::from_bits(0x3FBBE39EBE6F07C3, 0x3C5F7B8F29A05987), // 1.089419569898657998418608611626e-1
    F64x2::from_bits(0x3FBDDD21701EBA6E, 0x3C594EFFCD76FE58), // 1.166554354410693549660528495093e-1
    F64x2::from_bits(0x3FBFD5BA9AAC2F6E, 0xBC4CD37686760C17), // 1.243549945467614350313548491639e-1
    F64x2::from_bits(0x3FC0E6ADCCF40882, 0xBC6D71A31BB98D0D), // 1.320397616146387492746844065266e-1
    F64x2::from_bits(0x3FC1E1FAFB043727, 0xBC4B485914DACF8C), // 1.397088742891636451833677767391e-1
    F64x2::from_bits(0x3FC2DCBDB2FBA1FF, 0x3C58F28705561534), // 1.473614810886516356098027603968e-1
    F64x2::from_bits(0x3FC3D6EEE8C6626C, 0x3C661A3B0CE9281B), // 1.549967419239409823037143749335e-1
    F64x2::from_bits(0x3FC4D087A9DA4F17, 0x3C61F323F1ADF158), // 1.626138285979485753736415637616e-1
    F64x2::from_bits(0x3FC5C9811E3EC26A, 0xBC5054AB2C010F3D), // 1.702119252854744044904966070998e-1
    F64x2::from_bits(0x3FC6C1D4898933D9, 0xBC52954A7603C427), // 1.777902289926760707966247992158e-1
    F64x2::from_bits(0x3FC7B97B4BCE5B02, 0x3C5347B0B4F881CA), // 1.853479499956947648860259612285e-1
    F64x2::from_bits(0x3FC8B06EE2879C29, 0xBC6118CD30308C4F), // 1.928843122579746641970587106902e-1
    F64x2::from_bits(0x3FC9A6A8E96C8626, 0x3C4CF601E7B4348E), // 2.003985538258785146539457850344e-1
    F64x2::from_bits(0x3FCA9C231B403279, 0x3C60E8BBE89CCA85), // 2.078899272022629936053349831030e-1
    F64x2::from_bits(0x3FCB90D7529260A2, 0x3C217B10D2E0E5AB), // 2.153576996977380480244596271665e-1
    F64x2::from_bits(0x3FCC84BF8A742E6E, 0xBC595BDD0682EA26), // 2.228011537593945157710321221404e-1
    F64x2::from_bits(0x3FCD77D5DF205736, 0x3C6C648D1534597E), // 2.302195872768437302401709596798e-1
    F64x2::from_bits(0x3FCE6A148E96EC4D, 0x3C6866B22029F765), // 2.376123138654712524738836343256e-1
    F64x2::from_bits(0x3FCF5B75F92C80DD, 0x3C68AB6E3CF7AFBD), // 2.449786631268641541720824812113e-1
    F64x2::from_bits(0x3FD025FA510665B6, 0xBC7672DF6832FA48), // 2.523179808864271872059075075477e-1
    F64x2::from_bits(0x3FD09DC597D86362, 0x3C762E47390CB865), // 2.596296294082575310299464431840e-1
    F64x2::from_bits(0x3FD1151A362431CA, 0xBC74DC8DC9077B9F), // 2.669129875874004339552348564916e-1
    F64x2::from_bits(0x3FD18BF5A30BF178, 0x3C630CA4748B1BF9), // 2.741674511196587975993718983422e-1
    F64x2::from_bits(0x3FD2025567E47C96, 0xBC61832328F4290E), // 2.813924326491784512863373252053e-1
    F64x2::from_bits(0x3FD278372057EF46, 0xBC7077CDD36DFC81), // 2.885873618940773956236114199582e-1
    F64x2::from_bits(0x3FD2ED987A823CFE, 0x3C6B91258EA012CA), // 2.957516857504315485808979994627e-1
    F64x2::from_bits(0x3FD362773707EBCC, 0xBC6963A544B672D8), // 3.028848683749714055605560945056e-1
    F64x2::from_bits(0x3FD3D6D129271134, 0x3C7137CA41CC958A), // 3.099863912468834453178920488852e-1
    F64x2::from_bits(0x3FD44AA436C2AF0A, 0xBC75D5E43C55B3BA), // 3.170557532091470098090155766745e-1
    F64x2::from_bits(0x3FD4BDEE586890E7, 0xBC6E4DC77C22A757), // 3.240924704898717069873818070036e-1
    F64x2::from_bits(0x3FD530AD9951CD4A, 0xBC62566480884082), // 3.310960767041320949443387877569e-1
    F64x2::from_bits(0x3FD5A2E0175E0F4E, 0x3C713B7A8F82E457), // 3.380661228368254816598223422680e-1
    F64x2::from_bits(0x3FD614840309CFE2, 0xBC7A725715711F00), // 3.450021772071051088676812869001e-1
    F64x2::from_bits(0x3FD685979F5FA6FE, 0xBC4257814D1ADA9C), // 3.519038254149647861988092323316e-1
    F64x2::from_bits(0x3FD6F61941E4DEF1, 0xBC7C63AAE6F6E918), // 3.587706702705722203959200639265e-1
    F64x2::from_bits(0x3FD7660752817502, 0xBC4DD11791CC7600), // 3.656023317069668733128101132587e-1
    F64x2::from_bits(0x3FD7D5604B63B3F7, 0x3C769C885C2B249A), // 3.723984466767542219236550382837e-1
    F64x2::from_bits(0x3FD84422B8DF95D7, 0x3C7D76A0299B41B6), // 3.791586690334418341599883683469e-1
    F64x2::from_bits(0x3FD8B24D394A1B25, 0x3C7B6D0BA3748FA8), // 3.858826693980737758976954846072e-1
    F64x2::from_bits(0x3FD91FDE7CD0C662, 0x3C71074188054B53), // 3.925701350118285951655766597158e-1
    F64x2::from_bits(0x3FD98CD5454D6B18, 0x3C79E6C988FD0A77), // 3.992207695752525656147166961589e-1
    F64x2::from_bits(0x3FD9F93066168002, 0xBC7C827047C9439A), // 4.058342930748040952326027825139e-1
    F64x2::from_bits(0x3FDA64EEC3CC23FD, 0xBC724DEC1B50B7FF), // 4.124104415973873068997912896671e-1
    F64x2::from_bits(0x3FDAD00F5422058B, 0x3C7FC4C33891D2E8), // 4.189489671335528684572821440931e-1
    F64x2::from_bits(0x3FDB3A911DA65C6C, 0x3C7AE187B1CA5040), // 4.254496373700422895422636051808e-1
    F64x2::from_bits(0x3FDBA473378624A5, 0x3C7519A1B46E4AFF), // 4.319122354723482121002482623680e-1
    F64x2::from_bits(0x3FDC0DB4C94EC9F0, 0xBC7CC1CE70934C34), // 4.383365598579578054456160492148e-1
    F64x2::from_bits(0x3FDC76550AAD71F9, 0xBC774B8BFF7043E4), // 4.447224239609393412483431930515e-1
    F64x2::from_bits(0x3FDCDE53432C1351, 0xBC7A2CFA4418F1AD), // 4.510696559885234763756392572822e-1
    F64x2::from_bits(0x3FDD45AEC9EC862B, 0x3C689421163EF92D), // 4.573780986703208202304474687767e-1
    F64x2::from_bits(0x3FDDAC670561BB4F, 0x3C7A2B7F222F65E2), // 4.636476090008061162142562314612e-1
    F64x2::from_bits(0x3FDE127B6B0744B0, 0xBC52B0986398D4AB), // 4.698780579756869317238837674290e-1
    F64x2::from_bits(0x3FDE77EB7F175A34, 0x3C70E53DC1BF3435), // 4.760693303227612340751004202615e-1
    F64x2::from_bits(0x3FDEDCB6D43F8435, 0xBC5FC976330884E4), // 4.822213242278537357239469689631e-1
    F64x2::from_bits(0x3FDF40DD0B541418, 0xBC6A3992DC382A23), // 4.883339510564055238671649607471e-1
    F64x2::from_bits(0x3FDFA45DD3029259, 0xBC7CA563DC28D8B5), // 4.944071350712753472270372523319e-1
    F64x2::from_bits(0x3FE0039C73C1A40C, 0xBC8B32C949C9D593), // 5.004408131472941140300005149792e-1
    F64x2::from_bits(0x3FE034B709250488, 0x3C78F9B38D855410), // 5.064349344830967542119594467475e-1
    F64x2::from_bits(0x3FE0657E94DB30D0, 0xBC7D5B495F6349E6), // 5.123894603107377066666010205843e-1
    F64x2::from_bits(0x3FE095F30861A590, 0xBC7121B20A15A9F3), // 5.183043636035779962067055995152e-1
    F64x2::from_bits(0x3FE0C6145B5B43DA, 0x3C5974FA13B5404F), // 5.241796287829132483216496175045e-1
    F64x2::from_bits(0x3FE0F5E28B67E295, 0x3BE311B17EC990D0), // 5.300152514237931323005627690087e-1
    F64x2::from_bits(0x3FE1255D9BFBD2A9, 0xBC52BDAEE1C0EE35), // 5.358112379604637002690850687077e-1
    F64x2::from_bits(0x3FE154859637646A, 0xBC84BA7C548BF3C3), // 5.415676053918449723975989999416e-1
    F64x2::from_bits(0x3FE1835A88BE7C13, 0x3C8C621CEC00C301), // 5.472843809874369739852207703128e-1
    F64x2::from_bits(0x3FE1B1DC87904285, 0xBC621E8C8AEF8F29), // 5.529616019940283210522390754782e-1
    F64x2::from_bits(0x3FE1E00BABDEFEB4, 0xBC5928DF287A668F), // 5.585993153435624359715082164017e-1
    F64x2::from_bits(0x3FE20DE813E823B2, 0xBC8791D753EBB744), // 5.641975773624976077917430252673e-1
    F64x2::from_bits(0x3FE23B71E2CC9E6A, 0x3C6C421C9F38224E), // 5.697564534829784433238348916656e-1
    F64x2::from_bits(0x3FE268A940696DA6, 0x3C5D1348A04C73CC), // 5.752760179561178307303390604992e-1
    F64x2::from_bits(0x3FE2958E59308E31, 0xBC709E73B0C6C087), // 5.807563535676703992032744750015e-1
    F64x2::from_bits(0x3FE2C2215E024466, 0xBC44B810DA3A4BE1), // 5.861975513563606443330871822473e-1
    F64x2::from_bits(0x3FE2EE628406CBCA, 0x3C8C5D5E9FF0CF8D), // 5.915997103351114331458526589590e-1
    F64x2::from_bits(0x3FE31A52048874BE, 0x3C840CAB87A7AC24), // 5.969629372154015360144938318679e-1
    F64x2::from_bits(0x3FE345F01CCE37BB, 0x3C81021137C71102), // 6.022873461349641816821226942042e-1
    F64x2::from_bits(0x3FE3713D0DF6C504, 0xBC54F789E031606D), // 6.075730583890224200170815974795e-1
    F64x2::from_bits(0x3FE39C391CD4171A, 0xBC82304331D8BF46), // 6.128202021652413251433846354957e-1
    F64x2::from_bits(0x3FE3C6E491C78DC5, 0xBC8E145094FD0BA7), // 6.180289122825617964551808283643e-1
    F64x2::from_bits(0x3FE3F13FB89E96F4, 0x3C7ECF8B492644F0), // 6.231993299340659309924753490604e-1
    F64x2::from_bits(0x3FE41B4AE06FEA41, 0x3C63D60A53277652), // 6.283316024340096592176448127739e-1
    F64x2::from_bits(0x3FE445065B795B56, 0xBC7F76D0163F79C8), // 6.334258829691445662686954830593e-1
    F64x2::from_bits(0x3FE46E727EFE4716, 0xBC639B9B1B844CC9), // 6.384823303544375687098458370249e-1
    F64x2::from_bits(0x3FE4978FA3269EE1, 0x3C72419A87F2A458), // 6.435011087932843868028092287173e-1
    F64x2::from_bits(0x3FE4C05E22DE94E5, 0xBC8C0AC1F09F2EDF), // 6.484823876423005464090698119163e-1
    F64x2::from_bits(0x3FE4E8DE5BB6EC04, 0x3C84A33DBEB3796C), // 6.534263411807619628638934113116e-1
    F64x2::from_bits(0x3FE51110ADC5ED81, 0x3C723DCD6832A63E), // 6.583331483847559997836305857186e-1
    F64x2::from_bits(0x3FE538F57B89061F, 0xBC81BB74ABDA520C), // 6.632029927060932553632543102383e-1
    F64x2::from_bits(0x3FE5608D29C70C34, 0x3C89939CF0DE8088), // 6.680360618560202023579847712093e-1
    F64x2::from_bits(0x3FE587D81F732FBB, 0xBC75E5C9D8C5A950), // 6.728325475937631893114013292615e-1
    F64x2::from_bits(0x3FE5AED6C5909517, 0x3C87312F714A9436), // 6.775926455199251921409679777758e-1
    F64x2::from_bits(0x3FE5D58987169B18, 0x3C60028E4BC5E7CA), // 6.823165548747480782564299817112e-1
    F64x2::from_bits(0x3FE5FBF0D0D5CC4A, 0xBC5B4CFD000B7158), // 6.870044783412450003144135183117e-1
    F64x2::from_bits(0x3FE6220D115D7B8E, 0xBC62B785350EE8C1), // 6.916566218531998629800663181104e-1
    F64x2::from_bits(0x3FE647DEB8E20B90, 0xBC5ECA04023A51CF), // 6.962731944080235923905100965444e-1
    F64x2::from_bits(0x3FE66D663923E087, 0xBC76EA6FEBE8BBBA), // 7.008544078844501724579512817868e-1
    F64x2::from_bits(0x3FE692A40556FB6A, 0x3C8D94B95A8EA2CC), // 7.054004768650490822204263019559e-1
    F64x2::from_bits(0x3FE6B798920B3D99, 0xBC8A80386188C50E), // 7.099116184635248611916111509362e-1
    F64x2::from_bits(0x3FE6DC44551553AF, 0xBC5BF8863573828E), // 7.143880521567690317496541110864e-1
    F64x2::from_bits(0x3FE700A7C5784634, 0xBC78C34D25AADEF6), // 7.188299996216245054170141515259e-1
    F64x2::from_bits(0x3FE724C35B4FAE7B, 0x3C5948B32DB3499B), // 7.232376845763178795792637634429e-1
    F64x2::from_bits(0x3FE748978FBA8E0F, 0x3C47B2A6165884A1), // 7.276113326265106787829526909500e-1
    F64x2::from_bits(0x3FE76C24DCC6C6C0, 0x3C81952551ADC83D), // 7.319511711159165961691081923070e-1
    F64x2::from_bits(0x3FE78F6BBD5D315E, 0x3C8406A089803740), // 7.362574289814281317428352710891e-1
    F64x2::from_bits(0x3FE7B26CAD2E50FE, 0xBC8CE80DF30411FB), // 7.405303366126926909825420855292e-1
    F64x2::from_bits(0x3FE7D528289FA093, 0x3C8560821E2F3AA9), // 7.447701257160751857639310909741e-1
    F64x2::from_bits(0x3FE7F79EACB97898, 0x3C8FD5CA80EAD221), // 7.489770291829414558452980730009e-1
    F64x2::from_bits(0x3FE819D0B7158A4D, 0xBC7BF76229D3B917), // 7.531512809621943895247393702690e-1
    F64x2::from_bits(0x3FE83BBEC5CDEE22, 0x3C63107104FFC6C3), // 7.572931159369924530266468114325e-1
    F64x2::from_bits(0x3FE85D69576CC2C5, 0x3C66B66E7FC8B8C3), // 7.614027698055784264231855420876e-1
    F64x2::from_bits(0x3FE87ED0EADC5A2A, 0x3C70AF5AD957F4BC), // 7.654804789661444764216520449312e-1
    F64x2::from_bits(0x3FE89FF5FF57F1F8, 0xBC855B9A5E177A1B), // 7.695264804056582604068200359857e-1
    F64x2::from_bits(0x3FE8C0D9145CF49D, 0x3C8BEA4076DC4333), // 7.735410115925735392764068852830e-1
    F64x2::from_bits(0x3FE8E17AA99CC05E, 0xBC7EC182AB042F61), // 7.775243103733477667249308161224e-1
    F64x2::from_bits(0x3FE901DB3EEEF187, 0x3C868665E5603C8F), // 7.814766148726883079281840634961e-1
    F64x2::from_bits(0x3FE921FB54442D18, 0x3C81A62633145C07), // 7.853981633974483096156608458199e-1
];

// GENERATE: atanh_table F64x2 7 120
// ATANH_TBL[i] = atanh(i / 128)
static ATANH_TBL: [F64x2; 121] = [
    F64x2::from_bits(0x0000000000000000, 0x0000000000000000), // 0
    F64x2::from_bits(0x3F8000155588891B, 0xBC1645E17581DD13), // 7.812658951540420910323471276040e-3
    F64x2::from_bits(0x3F9000555888AD1D, 0xBC39C7CB16A754C5), // 1.562627175205221137920177875164e-2
    F64x2::from_bits(0x3F980120184F3DED, 0xBC396797D97ECFBB), // 2.344179294942521017357300574548e-2
    F64x2::from_bits(0x3FA001558891AEE2, 0x3C42D2774FF716DA), // 3.126017849066699476401224517265e-2
    F64x2::from_bits(0x3FA4029B471650A5, 0xBC494E8C5B9FB29C), // 3.908238642466813492618578268555e-2
    F64x2::from_bits(0x3FA8048185694820, 0xBC4096674E2A0E26), // 4.690937760882740545934144646905e-2
    F64x2::from_bits(0x3FAC07289F7B9245, 0x3C325A4A2CF58DE0), // 5.474211642084418033130394403175e-2
    F64x2::from_bits(0x3FB005588AD375AD, 0xBC41A7676AD4E1CC), // 6.258157147700300712676502386221e-2
    F64x2::from_bits(0x3FB2079DC9754943, 0xBC517E9E43B01797), // 7.042871635829776498015138333178e-2
    F64x2::from_bits(0x3FB40A74799E283E, 0x3C5E082A84FE2B08), // 7.828453034576996579297384750651e-2
    F64x2::from_bits(0x3FB60DED24F86C90, 0x3C476536D5F4F044), // 8.614999916646781917240350173259e-2
    F64x2::from_bits(0x3FB81218741F5A6A, 0x3C42807E177679B7), // 9.402611575146982626142002215373e-2
    F64x2::from_bits(0x3FBA170731EF7B1E, 0x3C5AB436CA05944E), // 1.019138810074591064043418708292e-1
    F64x2::from_bits(0x3FBC1CCA4EE78E03, 0xBC42BB1103E5C13E), // 1.098143046033826349625909318513e-1
    F64x2::from_bits(0x3FBE2372E49BCE54, 0x3C5555C8ED40626C), // 1.177284057737833557666497629572e-1
    F64x2::from_bits(0x3FC015891C9EAEF7, 0x3C6A6519E83DAE46), // 1.256572141404530388425688652009e-1
    F64x2::from_bits(0x3FC119DCE19BDBAF, 0x3C6F3124691D0470), // 1.336017705541201330103452288171e-1
    F64x2::from_bits(0x3FC21EBD96730F38, 0xBC5A701ACE932CFE), // 1.415631279579601052147029731236e-1
    F64x2::from_bits(0x3FC324342C318E81, 0xBC6F1497852F88F0), // 1.495423522747963006143171357998e-1
    F64x2::from_bits(0x3FC42A49B1BA196B, 0xBC670EF2144030E0), // 1.575405233199476850911799801319e-1
    F64x2::from_bits(0x3FC5310755C9FD18, 0xBC43D5E8218D1439), // 1.655587357417764838440147290404e-1
    F64x2::from_bits(0x3FC63876690E9070, 0x3C41A392DA64700C), // 1.735980999920943032261473821729e-1
    F64x2::from_bits(0x3FC740A0604C5ADB, 0x3C68AD317CDB57F2), // 1.816597433287004788474394929463e-1
    F64x2::from_bits(0x3FC8498ED69936DC, 0xBC690C066FCA62DA), // 1.897448108524518619777699951613e-1
    F64x2::from_bits(0x3FC9534B8FAAD565, 0xBC58F077AA182E9E), // 1.978544665813998461357113883324e-1
    F64x2::from_bits(0x3FCA5DE07A3B1BC2, 0x3C5CFF03EF7EF324), // 2.059898945646790387612582333063e-1
    F64x2::from_bits(0x3FCB6957B283EC91, 0x3C6448ED286E3F3B), // 2.141523000389935848228626954278e-1
    F64x2::from_bits(0x3FCC75BB84D40517, 0xBC284B1BDDC289C7), // 2.223429106307228391236070060177e-1
    F64x2::from_bits(0x3FCD8316703EB314, 0x3BAF9235DE46A021), // 2.305629776068590652860055922966e-1
    F64x2::from_bits(0x3FCE917329684475, 0x3C62A6AED3A61656), // 2.388137771781974419811200272975e-1
    F64x2::from_bits(0x3FCFA0DC9D713200, 0xBC3AB1DF878429E3), // 2.470966118584243517114963771238e-1
    F64x2::from_bits(0x3FD058AEFA811452, 0xBC7636250F41F7EF), // 2.554128118829953416027570481518e-1
    F64x2::from_bits(0x3FD0E1814BBD9D56, 0x3C679B88A2D962C4), // 2.637637366919610806511594050697e-1
    F64x2::from_bits(0x3FD16AEB1724557B, 0x3C7548148D78D887), // 2.721507764811900893800151389590e-1
    F64x2::from_bits(0x3FD1F4F2542FF85B, 0x3C57A5A117A66FEA), // 2.805753538267531983702614023279e-1
    F64x2::from_bits(0x3FD27F9D1EF3E177, 0xBC7D68A6061F2111), // 2.890389253875790565300052706136e-1
    F64x2::from_bits(0x3FD30AF1BA0717B8, 0x3C703246D5B0E323), // 2.975429836918652844495539815527e-1
    F64x2::from_bits(0x3FD396F69087FD7C, 0xBC70C4E3B04E71CF), // 3.060890590131390830421112081101e-1
    F64x2::from_bits(0x3FD423B2383A6343, 0x3C6E72370DD8C24A), // 3.146787213423076344709676180375e-1
    F64x2::from_bits(0x3FD4B12B73C1DD95, 0xBC3CF874EAB86927), // 3.233135824625262262193272012003e-1
    F64x2::from_bits(0x3FD53F6934FA63F8, 0xBC681D8B9AD374C7), // 3.319952981342448767642288069231e-1
    F64x2::from_bits(0x3FD5CE729F71680A, 0xBC71DEF95DFEF53B), // 3.407255703983770156887257188767e-1
    F64x2::from_bits(0x3FD65E4F0B01C08E, 0x3C53B2435A1D59F9), // 3.495061500061716939746119773168e-1
    F64x2::from_bits(0x3FD6EF060694F581, 0xBC7F8FB0B9C0BFB9), // 3.583388389850697134861222664901e-1
    F64x2::from_bits(0x3FD7809F5B0CB028, 0xBC4E3C176344A2FF), // 3.672254933505905117290140825748e-1
    F64x2::from_bits(0x3FD813230E574D59, 0x3C7D6EE17087C874), // 3.761680259751379573558769481235e-1
    F64x2::from_bits(0x3FD8A69966B2D128, 0xBC799ED712A2C701), // 3.851684096255376443629452311027e-1
    F64x2::from_bits(0x3FD93B0AEE21C2C8, 0x3C41EB08151FBD93), // 3.942286801821350847305921223695e-1
    F64x2::from_bits(0x3FD9D0807615C643, 0xBC53140EEADF73FC), // 4.033509400534038239191750345804e-1
    F64x2::from_bits(0x3FDA67031B542059, 0xBC78CC9C3D9CD60C), // 4.125373618012467054345414310910e-1
    F64x2::from_bits(0x3FDAFE9C4A18B0E3, 0x3C58CDBD101F4B14), // 4.217901919935355735052540543065e-1
    F64x2::from_bits(0x3FDB9755C27C59DF, 0x3C7AE8C2BEEA8DFE), // 4.311117553019396486911111875917e-1
    F64x2::from_bits(0x3FDC31399D243E72, 0x3C699DFAD9C3F0EE), // 4.405044588647576531211615676291e-1
    F64x2::from_bits(0x3FDCCC52503FC6FE, 0xBC6DB36DCEE185BE), // 4.499707969363127983953722569885e-1
    F64x2::from_bits(0x3FDD68AAB4DBE74B, 0xBC7788B596EC6913), // 4.595133558465149055272627414548e-1
    F64x2::from_bits(0x3FDE064E0C92C396, 0xBC6CA9772FC5EAAF), // 4.691348192964651540081311789826e-1
    F64x2::from_bits(0x3FDEA548079F8314, 0x3C3B8ABD5472DB5B), // 4.788379740185046988197002358723e-1
    F64x2::from_bits(0x3FDF45A4CB5EE467, 0x3C7D159F0CC3EBFA), // 4.886257158319211288091487423978e-1
    F64x2::from_bits(0x3FDFE770F9460541, 0xBC6C1E34E2AF3E1C), // 4.985010561286636210547400635528e-1
    F64x2::from_bits(0x3FE0455CDB2CE279, 0x3C651A7E0FF217EB), // 5.084671288269212532857080259495e-1
    F64x2::from_bits(0x3FE097C659991EC2, 0xBC89FBCBB91EE41C), // 5.185271978343381598105082289235e-1
    F64x2::from_bits(0x3FE0EAFC1A4B81EA, 0x3C8206FC6402695A), // 5.286846650670302995762989751132e-1
    F64x2::from_bits(0x3FE13F058E611D13, 0xBC8E7E14EFA7917C), // 5.389430790754963810978515887054e-1
    F64x2::from_bits(0x3FE193EA7AAD030B, 0xBC8A256F99CAABEB), // 5.493061443340548456976226184613e-1
    F64x2::from_bits(0x3FE1E9B2FD18D91B, 0x3C80B8E4344FA6FB), // 5.597777312566764319795439180474e-1
    F64x2::from_bits(0x3FE240679272D92C, 0x3C70C7B60E19F82A), // 5.703618870091182880937786485808e-1
    F64x2::from_bits(0x3FE298111CB3F8AC, 0x3C8E2C8F7949405B), // 5.810628471952177766490514541442e-1
    F64x2::from_bits(0x3FE2F0B8E9CA2471, 0xBC866259A0133301), // 5.918850485042082721073504583688e-1
    F64x2::from_bits(0x3FE34A68BAF3E921, 0x3C7B7D83CD1CE65B), // 6.028331424161345280818734777216e-1
    F64x2::from_bits(0x3FE3A52ACCBC786F, 0x3C71BFFDB9B19E0C), // 6.139120100740579498345810097241e-1
    F64x2::from_bits(0x3FE40109DFA8CCB5, 0x3C6CAA0D2A74F2E1), // 6.251267784449711220029243115291e-1
    F64x2::from_bits(0x3FE45E1141A8C010, 0xBC8F9CDE95328C5E), // 6.364828379064437220480829615046e-1
    F64x2::from_bits(0x3FE4BC4CD8614BF5, 0x3C7FC876CBA989E9), // 6.479858614133024162013071610475e-1
    F64x2::from_bits(0x3FE51BC92C67DFA3, 0x3C79525165EA055B), // 6.596418254184651883492312780727e-1
    F64x2::from_bits(0x3FE57C937589DD3F, 0x3C898FD7AC1BAC64), // 6.714570327448327490720778183880e-1
    F64x2::from_bits(0x3FE5DEB9A83EE95C, 0x3C7C2DD299AAC058), // 6.834381376313945177956462066152e-1
    F64x2::from_bits(0x3FE6424A8468E3D1, 0xBC8CC43A5853BD18), // 6.955921732070412034113714257340e-1
    F64x2::from_bits(0x3FE6A755A5893549, 0xBC7E58685ABD484F), // 7.079265816807175107835418243005e-1
    F64x2::from_bits(0x3FE70DEB9498B947, 0xBC89F7D0D7D2A389), // 7.204492475773712316932689093844e-1
    F64x2::from_bits(0x3FE7761DDBB61597, 0x3C8968AE8146FDEA), // 7.331685343967135223291211023214e-1
    F64x2::from_bits(0x3FE7DFFF1BE5F384, 0xBC87E7CBA6A8D6D1), // 7.460933251273762086250101700404e-1
    F64x2::from_bits(0x3FE84BA3252982B5, 0x3C812B2596ABBB4D), // 7.592330671141868401473022408908e-1
    F64x2::from_bits(0x3FE8B91F113A34A3, 0xBC82510CBFF5976D), // 7.725978218528733691753532291915e-1
    F64x2::from_bits(0x3FE9288961460ABE, 0xBC49EDB6A2E96A7A), // 7.861983203768756450412627805320e-1
    F64x2::from_bits(0x3FE999FA1F16860C, 0x3C805B1AAC05D48E), // 8.000460250079313446259286280449e-1
    F64x2::from_bits(0x3FEA0D8B021DC00C, 0xBC72211531505AE2), // 8.141531983692415884019004183288e-1
    F64x2::from_bits(0x3FEA835798FA0CBA, 0x3C8017A03CDAEBF9), // 8.285329807116774971034192367673e-1
    F64x2::from_bits(0x3FEAFB7D78197BEC, 0x3C86D033540BA421), // 8.431994767851143497924881887707e-1
    F64x2::from_bits(0x3FEB761C6E44955F, 0xBC87765DAA3A32DD), // 8.581678537054065617928921947488e-1
    F64x2::from_bits(0x3FEBF356BFFBEDC0, 0x3C7DAE5B415BEC73), // 8.734544515313516201841445101607e-1
    F64x2::from_bits(0x3FEC73516AC03329, 0xBC550DDB1592ABED), // 8.890769085861381880596292116648e-1
    F64x2::from_bits(0x3FECF6347191F5B6, 0xBC8517748DEFFC16), // 9.050543039481257679361062186731e-1
    F64x2::from_bits(0x3FED7C2B3438300B, 0xBC723102F6ADC943), // 9.214073200141695682984549913845e-1
    F64x2::from_bits(0x3FEE0564D32D9391, 0x3C8A26ED45A8E606), // 9.381584286280590277569031450404e-1
    F64x2::from_bits(0x3FEE9214A278F736, 0xBC8B5C0AC0DB1011), // 9.553321049968192615894581974177e-1
    F64x2::from_bits(0x3FEF2272AE325A57, 0x3C851BDA525B3C98), // 9.729550745276566525526763717216e-1
    F64x2::from_bits(0x3FEFB6BC5412C9C8, 0x3C71510A70EE438D), // 9.910565988596369430314224078775e-1
    F64x2::from_bits(0x3FF0279A7B19BE98, 0x3C9BEE52F27876BB), // 1.009668808805065251293063461734e0
    F64x2::from_bits(0x3FF07613660E2537, 0xBC94345A063CB619), // 1.028827093747464323915054772236e0
    F64x2::from_bits(0x3FF0C6F1EC420786, 0x3C737DE6BEC7A18A), // 1.048570559389618423149459962684e0
    F64x2::from_bits(0x3FF11A6385E30CF5, 0x3C9B22FF473CAEE8), // 1.068942568774955287127749946583e0
    F64x2::from_bits(0x3FF1709ACE96EE39, 0xBC8FF8FE086ECC63), // 1.089991385450856664677012361736e0
    F64x2::from_bits(0x3FF1C9D0540158C9, 0x3C90AACEADA71A35), // 1.111770942826796395680511445939e0
    F64x2::from_bits(0x3FF226438E777C99, 0x3C9C6C5FF4A16C83), // 1.134341770659182167894013397720e0
    F64x2::from_bits(0x3FF2863C0EA8B8CF, 0xBC7F0006840A516A), // 1.157772118818275437247736021429e0
    F64x2::from_bits(0x3FF2EA0AEE5F5810, 0xBC80E5CEABFF59CE), // 1.182139330999692787391017168723e0
    F64x2::from_bits(0x3FF3520C9718089F, 0x3C791E94F0395701), // 1.207531538210367982460556011316e0
    F64x2::from_bits(0x3FF3BEAAF7978C58, 0xBC920D97E442CDF6), // 1.234049765735809538007645520247e0
    F64x2::from_bits(0x3FF430604AB1313F, 0xBC84444C8BB4666D), // 1.261810580984345362779581673258e0
    F64x2::from_bits(0x3FF4A7BA9E66A99B, 0xBC7E335D57E61AF8), // 1.290949457887655001282198443736e0
    F64x2::from_bits(0x3FF525605D6FC1B9, 0x3C995A2861C9EBE8), // 1.321625103937647293356814488194e0
    F64x2::from_bits(0x3FF5AA16394D481F, 0x3C3341C89935864A), // 1.354025100551105032998002285074e0
    F64x2::from_bits(0x3FF636C6FDA710DD, 0x3C9EBF2B0B238900), // 1.388373366194222464735203374281e0
    F64x2::from_bits(0x3FF6CC8E17E5416C, 0xBC8D5F7D928E94AB), // 1.424940198270713891509267206244e0
    F64x2::from_bits(0x3FF76CC6077EFB49, 0x3C9F00C98BC72935), // 1.464056042939505860461369371524e0
    F64x2::from_bits(0x3FF8191C98CE5694, 0x3C9DEB16A663B4E0), // 1.506130787752600778678071933252e0
    F64x2::from_bits(0x3FF8D3AFF2A9E9A5, 0x3C948C8211B23EAC), // 1.551681468873178220374760621074e0
    F64x2::from_bits(0x3FF99F3986EE6739, 0x3C9078C142854C77), // 1.601373221469158560330624638367e0
    F64x2::from_bits(0x3FFA7F4FB6893C73, 0xBC8371A499649CE3), // 1.656081879645878906636012199804e0
    F64x2::from_bits(0x3FFB78CE48912B5A, 0xBC6DD1A2BB8F2588), // 1.716993602242573122964582162271e0
];

// GENERATE: bessel_debye_poly 14
// P(a, b) = sum(C[k][i] * a^(k - i) * b^i), for 1 <= k <= 14 and 0 <= i <= k
const DEBYE_A1: F64x2 = F64x2::from_bits(0xBFCAAAAAAAAAAAAB, 0x3C65555555555555); // -2.083333333333333333333333333333e-1
const DEBYE_B1: F64x2 = F64x2::from_bits(0x3FC0000000000000, 0x0000000000000000); // 1.250000000000000000000000000000e-1
// DEBYE_POLY has C[m + i][i] for a^m * b^i, with m + i >= 2, grouped by i
static DEBYE_POLY: [f64; 120] = [
    // b^0
    f64::from_bits(0x0000000000000000), // 0e0
    f64::from_bits(0x0000000000000000), // 0e0
    f64::from_bits(0x3FD4000000000000), // 3.125e-1
    f64::from_bits(0xBFEEB1C71C71C71C), // -9.592013888888888e-1
    f64::from_bits(0x4011A80000000000), // 4.4140625e0
    f64::from_bits(0xC03AF61555555555), // -2.6961263020833332e1
    f64::from_bits(0x40699E5555555555), // 2.0494791666666666e2
    f64::from_bits(0xC09D1C4901E79E7A), // -1.86307129632859e3
    f64::from_bits(0x40D33F3DD4000000), // 1.9708966064453125e4
    f64::from_bits(0xC10D08D59B7F097B), // -2.3785070092589766e5
    f64::from_bits(0x41489AF1EC280000), // 3.225059844970703e6
    f64::from_bits(0xC18725A9E17C80FD), // -4.85430361857929e7
    f64::from_bits(0x41C7EFD415720955), // 8.031867308909098e8
    f64::from_bits(0xC20AFD86F1DAD055), // -1.4490328635351725e10
    f64::from_bits(0x42507A81A9453526), // 2.8309917826083044e11
    // b^1
    f64::from_bits(0x0000000000000000), // 0e0
    f64::from_bits(0xBFD8000000000000), // -3.75e-1
    f64::from_bits(0x3FFBA00000000000), // 1.7265625e0
    f64::from_bits(0xC025300000000000), // -1.059375e1
    f64::from_bits(0x4054389000000000), // 8.08837890625e1
    f64::from_bits(0xC0870E8000000000), // -7.378125e2
    f64::from_bits(0x40BE90E642000000), // 7.824899444580078e3
    f64::from_bits(0xC0F718B098000000), // -9.4603037109375e4
    f64::from_bits(0x41339929C8F5C000), // 1.2843937849998474e6
    f64::from_bits(0xC1727435711E0000), // -1.935035906982422e7
    f64::from_bits(0x41B318AC26D38404), // 3.2038403882623315e8
    f64::from_bits(0xC1F58B0BACE6A200), // -5.782944462414551e9
    f64::from_bits(0x423A50C9F89BBE53), // 1.1302456335574345e11
    f64::from_bits(0xC2814D6E8B5577CF), // -2.378033097390976e12
    // b^2
    f64::from_bits(0x3FB0000000000000), // 6.25e-2
    f64::from_bits(0xBFEA8CCCCCCCCCCD), // -8.296875e-1
    f64::from_bits(0x4020980000000000), // 8.296875e0
    f64::from_bits(0xC055EC6000000000), // -8.7693359375e1
    f64::from_bits(0x408FED2000000000), // 1.021640625e3
    f64::from_bits(0xC0C9C06239000000), // -1.3184767364501953e4
    f64::from_bits(0x4106ECCF36000000), // 1.878019013671875e5
    f64::from_bits(0xC14664CF26968000), // -2.9351983014678955e6
    f64::from_bits(0x4187DAE6F42A8000), // 5.002774252075195e7
    f64::from_bits(0xC1CB8CF4BFABD69E), // -9.244450553424871e8
    f64::from_bits(0x42112802A6EB9B00), // 1.8421426618901367e10
    f64::from_bits(0xC256EE4A13BDE214), // -3.939484505515325e11
    f64::from_bits(0x42A05FF2CF31A19E), // 9.002140801232809e12
    // b^3
    f64::from_bits(0x3FB0AAAAAAAAAAAB), // 6.510416666666667e-2
    f64::from_bits(0xC001C00000000000), // -2.21875e0
    f64::from_bits(0x40444578E38E38E4), // 4.0542751736111114e1
    f64::from_bits(0xC0850D0000000000), // -6.73625e2
    f64::from_bits(0x40C614546A99999A), // 1.1304659503173829e4
    f64::from_bits(0xC108358E74000000), // -1.98321806640625e5
    f64::from_bits(0x414C1A8C91D0D555), // 3.6836091391855874e6
    f64::from_bits(0xC191583D58060000), // -7.274888600585938e7
    f64::from_bits(0x41D6C71E6F696DFC), // 1.528592829647338e9
    f64::from_bits(0xC21FCA603AC73000), // -3.4134822577796875e10
    f64::from_bits(0x4267888B002C8A67), // 8.086004698603251e11
    f64::from_bits(0xC2B27082C6510086), // -2.0274439672064523e13
    // b^4
    f64::from_bits(0x3FBA000000000000), // 1.015625e-1
    f64::from_bits(0xC01BEDDB6DB6DB6E), // -6.982282366071429e0
    f64::from_bits(0x406A340000000000), // 2.09625e2
    f64::from_bits(0xC0B4143616000000), // -5.140211273193359e3
    f64::from_bits(0x40FD11FF3E000000), // 1.1907195263671875e5
    f64::from_bits(0xC14500B29B0FC000), // -2.752869211418152e6
    f64::from_bits(0x418F0C429E5D0000), // 6.5112147795410156e7
    f64::from_bits(0xC1D7C334CE2B51A8), // -1.594676024676859e9
    f64::from_bits(0x4222F1A0C97E7750), // 4.068162886323303e10
    f64::from_bits(0xC26F8C944B7E2182), // -1.0840201205610471e12
    f64::from_bits(0x42BB77CD6CBF55CC), // 3.0201361514325797e13
    // b^5
    f64::from_bits(0x3FCAD33333333333), // 2.095703125e-1
    f64::from_bits(0xC039500000000000), // -2.53125e1
    f64::from_bits(0x4092235158000000), // 1.1608294372558594e3
    f64::from_bits(0xC0E3934E50000000), // -4.0090447265625e4
    f64::from_bits(0x4132EC3CE1D48000), // 1.2401248821487427e6
    f64::from_bits(0xC1819295F64B999A), // -3.6852414786914065e7
    f64::from_bits(0x41D0382F56529ADA), // 1.0884703612907014e9
    f64::from_bits(0xC21E4E387E0AFF00), // -3.2540336002749023e10
    f64::from_bits(0x426CF3F1FF1223FB), // 9.948144580011244e11
    f64::from_bits(0xC2BC7461C55BCF71), // -3.128618210811144e13
    // b^6
    f64::from_bits(0x3FE12AAAAAAAAAAB), // 5.364583333333334e-1
    f64::from_bits(0xC059FE7580000000), // -1.0397592163085938e2
    f64::from_bits(0x40BAF846C0000000), // 6.9042763671875e3
    f64::from_bits(0xC113C48E43595555), // -3.238755657704671e5
    f64::from_bits(0x4168D88CB9840000), // 1.3026405797363281e7
    f64::from_bits(0xC1BD0D79CC455C56), // -4.874224122709402e8
    f64::from_bits(0x4210782A608CB700), // 1.768387997117871e10
    f64::from_bits(0xC262878855ED1D02), // -6.366661630489065e11
    f64::from_bits(0x42B4F9D4C6045D6C), // 2.306324915926142e13
    // b^7
    f64::from_bits(0x3FFA358492492492), // 1.6380658830915178e0
    f64::from_bits(0xC07DD9D800000000), // -4.77615234375e2
    f64::from_bits(0x40E584B925600000), // 4.406978581237793e4
    f64::from_bits(0xC144CF134BC00000), // -2.727462591796875e6
    f64::from_bits(0x41A0A246965850B0), // 1.395351791724906e8
    f64::from_bits(0xC1F811F4FB3A1800), // -6.461280179630859e9
    f64::from_bits(0x425086BF0482528C), // 2.839205688412898e11
    f64::from_bits(0xC2A621A72EBE299E), // -1.2166897295124809e13
    // b^8
    f64::from_bits(0x40173D4000000000), // 5.809814453125e0
    f64::from_bits(0xC0A2F693BAAE8BA3), // -2.4272885336442428e3
    f64::from_bits(0x4112620C18400000), // 3.011870236816406e5
    f64::from_bits(0xC176E64DA48BEA40), // -2.4011994284158945e7
    f64::from_bits(0x41D6DB88D3CE8E00), // 1.533944655227417e9
    f64::from_bits(0xC2342A5623638572), // -8.660963414752127e10
    f64::from_bits(0x429090DF9A529618), // 4.553603191973523e12
    // b^9
    f64::from_bits(0x403779A1F8E38E39), // 2.3475127749972874e1
    f64::from_bits(0xC0CA6B7F70000000), // -1.352699560546875e4
    f64::from_bits(0x4140C43C6016039A), // 2.197624750671816e6
    f64::from_bits(0xC1AA5EFAEC20EAAB), // -2.2121611806429037e8
    f64::from_bits(0x42102D6797168F04), // 1.7370310085639664e10
    f64::from_bits(0xC2712E3C85439FB2), // -1.180642595897981e12
    // b^10
    f64::from_bits(0x405A9DF0CCCCCCCD), // 1.06467822265625e2
    f64::from_bits(0xC0F407E90D82E762), // -8.204656579866775e4
    f64::from_bits(0x417046EB88980000), // 1.7067704537109375e7
    f64::from_bits(0xC1DFC7170A6881B8), // -2.1325650336329174e9
    f64::from_bits(0x4247A318C40475EE), // 2.0304111412092133e11
    // b^11
    f64::from_bits(0x4080BD1FC8B1745D), // 5.35640519510616e2
    f64::from_bits(0xC1206C681B000000), // -5.38164052734375e5
    f64::from_bits(0x41A0C4EAA0E4E02E), // 1.4067028844702286e8
    f64::from_bits(0xC21405A4B5FD9DC0), // -2.1498506623404053e10
    // b^12
    f64::from_bits(0x40A7217AC5555555), // 2.960739786783854e3
    f64::from_bits(0xC14CF6A1EE88FDC0), // -3.79629186355564e6
    f64::from_bits(0x41D2481AA6D85100), // 1.2268612113799438e9
    // b^13
    f64::from_bits(0x40D16B51E66C789E), // 1.7837279688947478e4
    f64::from_bits(0xC17B55708CB18000), // -2.866151279333496e7
    // b^14
    f64::from_bits(0x40FC6733F2524925), // 1.1633924666050503e5
];

// `DEBYE_POLY` has the coefficients of the powers of `b` up to `DEBYE_TERMS`
// (see `debye_poly`)
const _: () = assert!(DEBYE_POLY.len() == (DEBYE_TERMS + 1) * (DEBYE_TERMS + 2) / 2);

/// Returns `true` if `R^3 / n^2 >= DEBYE_DIST_MIN`, where
/// `R = sqrt(|n^2 - x^2|)`, so the Debye expansions can be used for the
/// order `n` and the argument `x` (see `DEBYE_DIST_MIN`).
#[inline]
fn debye_valid(n: f64, x: f64) -> bool {
    // The rounding errors only move the limit slightly
    let r2 = ((n - x) * (n + x)).abs();
    (r2 * r2) * r2 >= (DEBYE_DIST_MIN * DEBYE_DIST_MIN) * ((n * n) * (n * n))
}

/// Returns the terms of degree 2 to `DEBYE_TERMS` of
/// `P(a, b) = sum(C[k][i] * a^(k - i) * b^i)` (see `DEBYE_POLY`), with plain
/// `f64` arithmetic.
///
/// The terms of degree `k` are `W[k](p) / n^k` with `a = p^3 / n` and
/// `b = p / n`, where `ln(sum(U[k](p) / n^k)) = sum(W[k](p) / n^k)` and `U[k]`
/// are the polynomials of the Debye expansions.
#[inline]
fn debye_poly(a: f64, b: f64) -> f64 {
    // sum(R[i](a) * b^i), where R[i] has the coefficients of `b^i`
    let mut acc = 0.0;
    let mut end = DEBYE_POLY.len();
    for i in (0..=DEBYE_TERMS).rev() {
        let start = end - (DEBYE_TERMS + 1 - i);
        let r = DEBYE_POLY[start..end]
            .iter()
            .rev()
            .fold(0.0, |s, &c| s * a + c);
        acc = acc * b + r;
        end = start;
    }
    acc
}

/// Returns `(re, im)` such that `re + i * im` is the sum of the terms of
/// degree 2 to `DEBYE_TERMS` of `P(i * a, -i * b)` (see `debye_poly`), with
/// plain `f64` arithmetic.
#[inline]
fn debye_poly_complex(a: f64, b: f64) -> (f64, f64) {
    // With `z = re + i * im`, z * (i * a) = -im * a + i * re * a and
    // z * (-i * b) = im * b - i * re * b
    let (mut acc_re, mut acc_im) = (0.0, 0.0);
    let mut end = DEBYE_POLY.len();
    for i in (0..=DEBYE_TERMS).rev() {
        let start = end - (DEBYE_TERMS + 1 - i);
        let (r_re, r_im) = DEBYE_POLY[start..end]
            .iter()
            .rev()
            .fold((0.0, 0.0), |(re, im), &c| (c - im * a, re * a));
        (acc_re, acc_im) = (r_re + acc_im * b, r_im - acc_re * b);
        end = start;
    }
    (acc_re, acc_im)
}

/// Returns `u * (1 + s / 3 + s^2 / 5 + ... + s^9 / 19)`, with `s = u^2`
/// (`hyperbolic = true`) or `s = -u^2`, which is `atanh(u)` or `atan(u)` for
/// `|u| <= 2^-4.9`, with a relative error of about 2^-100 with
/// `|u| <= 2^-6.8`, and 2^-85 with larger ones.
///
/// The first terms are calculated as sums of two `f64`, and the ones from
/// `s^3 / 7` (less than 2^-43.5 with `|u| <= 2^-6.8`, and 2^-32.2 with
/// larger ones) with plain `f64` arithmetic. The terms that are left out are
/// less than 2^-102 relative to the result.
#[inline]
fn odd_series(u: F64x2, hyperbolic: bool) -> F64x2 {
    let s = u * u;
    let s = if hyperbolic { s } else { -s };
    u * (1.0 + s * odd_series_tail(s))
}

/// Returns `1/3 + s / 5 + s^2 / 7 + ... + s^8 / 19` as a sum of two `f64`,
/// which is the factor of `u^3` in `odd_series` (see there).
#[inline]
fn odd_series_tail(s: F64x2) -> F64x2 {
    let sh = s.hi();
    let tail_hi = (1.0 / 13.0) + sh * ((1.0 / 15.0) + sh * ((1.0 / 17.0) + sh * (1.0 / 19.0)));
    let tail = (1.0 / 7.0) + sh * ((1.0 / 9.0) + sh * ((1.0 / 11.0) + sh * tail_hi));
    let q = FRAC_1_5 + sh * tail;
    FRAC_1_3 + s * q
}

/// Returns `(i, i / 128)`, where `i / 128` is the closest to `y`, for
/// `0 <= y <= 1`.
#[inline]
fn table_point(y: f64) -> (usize, f64) {
    let i = (y * 128.0 + 0.5) as usize;
    (i, i as f64 / 128.0)
}

/// Returns `atan(y)` as a sum of two `f64`, for `0 <= y <= 1`, with a
/// relative error of about 2^-100.
///
/// `atan(y) = atan(c) + atan(u)`, where `c = i / 128` is the closest to `y`
/// (see `ATAN_TBL`) and `u = (y - c) / (1 + y * c)`, so `|u| <= 2^-8`.
#[inline]
fn atan_precise(y: F64x2) -> F64x2 {
    let (i, c) = table_point(y.hi());
    // `y.hi() - c` is exact (by Sterbenz lemma, or because `c = 0`)
    let d = F64x2::add11(y.hi() - c, y.lo());
    let p = F64x2::mul11(y.hi(), c) + y.lo() * c;
    ATAN_TBL[i] + odd_series(d / (1.0 + p), false)
}

/// Returns `atanh(w)` as a sum of two `f64`, for `0 <= w <= 0.9375`, with a
/// relative error of about 2^-100 with `w <= 0.75`, and 2^-90 with larger
/// ones.
///
/// `atanh(w) = atanh(c) + atanh(u)`, where `c = i / 128` is the closest to
/// `w` (see `ATANH_TBL`) and `u = (w - c) / (1 - w * c)`, so
/// `|u| <= 2^-8 / (1 - 0.75^2) < 2^-6.8` with `w <= 0.75`, and
/// `|u| <= 2^-8 / (1 - 0.9375^2) < 2^-4.9` otherwise. Then the error of
/// `atanh(u)` (2^-85 relative to it, see `odd_series`) is less than 2^-89.9,
/// while `atanh(w) > 0.97`.
#[inline]
fn atanh_precise(w: F64x2) -> F64x2 {
    let (i, c) = table_point(w.hi());
    // `w.hi() - c` is exact (by Sterbenz lemma, or because `c = 0`)
    let d = F64x2::add11(w.hi() - c, w.lo());
    let p = F64x2::mul11(w.hi(), c) + w.lo() * c;
    ATANH_TBL[i] + odd_series(d / (1.0 - p), true)
}

/// Returns `R = sqrt(|n^2 - x^2|)` as a sum of two `f64`, with a relative
/// error of about 2^-100 (see `sqrt_parts`), for `x != n`.
#[inline]
fn debye_r(n: f64, x: f64) -> F64x2 {
    // |n^2 - x^2| = |n - x| * (n + x), where both factors are exact sums of
    // two `f64`
    let d = F64x2::sub11(n, x).abs();
    let r2 = d * F64x2::add11(n, x);
    let (hi, lo) = sqrt_parts(r2.hi(), r2.lo());
    F64x2::fast_add11(hi, lo)
}

/// Returns `(atanh(y) - y) / y^3` (`hyperbolic = true`) or
/// `(y - atan(y)) / y^3` as a sum of two `f64`, for `0 <= y <= 0.75`, with a
/// relative error of about 2^-82.
///
/// When `y < 1/256`, it is calculated from the series of `odd_series` (see
/// `odd_series_tail`), with an error of about 2^-88. Otherwise, the
/// subtraction cancels all but a fraction of about `y^2 / 3 > 2^-17.6`, which
/// increases the relative error of `atanh_precise` and `atan_precise`
/// (2^-100) to about 2^-82.
#[inline]
fn odd_excess(y: F64x2, hyperbolic: bool) -> F64x2 {
    let s = y * y;
    if y.hi() < 1.0 / 256.0 {
        odd_series_tail(if hyperbolic { s } else { -s })
    } else if hyperbolic {
        (atanh_precise(y) - y) / (s * y)
    } else {
        (y - atan_precise(y)) / (s * y)
    }
}

/// Returns `cbrt(a)` as a sum of two `f64`, for `a.hi() >= 1`, with a
/// relative error of about 2^-103.
#[inline]
fn cbrt_parts(a: F64x2) -> F64x2 {
    // A Newton iteration from `y = cbrt(a.hi())` (with a relative error of
    // about 2^-53) squares its error:
    // cbrt(a) ~= y + (a - y^3) / (3 * y^2)
    // where `y^2` is exact and `a - y^3` has an error of about 2^-104 * a
    let y = crate::generic::cbrt(a.hi());
    let y2 = F64x2::mul11(y, y);
    let res = a - y2 * y;
    F64x2::fast_add11(y, res.hi() / (3.0 * y2.hi()))
}

/// Number of terms of the Taylor series of `airy`.
const AIRY_TERMS: usize = 22;

/// Number of terms of the Taylor series of `airy` that are calculated as sums
/// of two `f64`.
const AIRY_DW_TERMS: usize = 10;

/// Returns `(y, dy)` such that `y` is `Ai(t)` (`second = false`) or `Bi(t)`
/// as a sum of two `f64`, and `dy ~= Ai'(t)` or `Bi'(t)`, for the arguments
/// up to 1/8 beyond the first and the last points of `AIRY_TBL`.
///
/// They are calculated from the Taylor series at the closest point `t0` of
/// `AIRY_TBL`, of `h = t - t0` (`|h| <= 1/8`), whose coefficients follow
/// from `y'' = t * y`:
///
/// `(k + 2) * (k + 1) * c[k + 2] = t0 * c[k] + c[k - 1]`
///
/// Relative to the smallest magnitude of the function in the interval
/// (`t0 > 0`), or to its modulus `sqrt(Ai(t)^2 + Bi(t)^2)` (`t0 <= 0`), the
/// terms `c[k] * h^k` are less than 2^-25.9 from `k = AIRY_DW_TERMS`, which
/// are calculated with plain `f64` arithmetic, and the first one that is left
/// out is less than 2^-79.6. So the error of `y` is about 2^-77 relative to
/// the same values.
///
/// `dy` is calculated with plain `f64` arithmetic (its error is about 2^-51
/// relative to it or to the modulus of the derivatives), which is enough for
/// the term of `olver` where it is used.
#[inline]
fn airy(t: F64x2, second: bool) -> (F64x2, f64) {
    // The closest point, which can only be beyond the last one (by rounding
    // to the next index) when `t` is exactly 1/8 beyond it, or very close
    let i = (((t.hi() - AIRY_T_MIN) * 4.0 + 0.5) as usize).min(AIRY_TBL.len() - 1);
    let t0 = AIRY_T_MIN + 0.25 * i as f64;
    // `t.hi() - t0` is exact: `t0` is a multiple of 1/4, so of the ULP of
    // `t.hi()`, and the difference is at most `|t.hi()|` (it is `t.hi()` when
    // `t0 = 0`)
    let h = F64x2::add11(t.hi() - t0, t.lo());
    let row = &AIRY_TBL[i];
    let j = if second { 2 } else { 0 };

    let mut c = [F64x2::new1(0.0); AIRY_DW_TERMS];
    c[0] = row[j];
    c[1] = row[j + 1];
    c[2] = (c[0] * t0) * 0.5;
    for k in 3..AIRY_DW_TERMS {
        let s = c[k - 2] * t0 + c[k - 3];
        let m = (k * (k - 1)) as f64;
        c[k] = F64x2::div11(s.hi(), m) + s.lo() / m;
    }
    let mut f = [0.0; AIRY_TERMS];
    for k in 0..AIRY_DW_TERMS {
        f[k] = c[k].hi();
    }
    for k in AIRY_DW_TERMS..AIRY_TERMS {
        f[k] = (t0 * f[k - 2] + f[k - 3]) / (k * (k - 1)) as f64;
    }

    let hh = h.hi();
    let tail = f[AIRY_DW_TERMS..]
        .iter()
        .rev()
        .fold(0.0, |s, &ck| s * hh + ck);
    let y = c.iter().rev().fold(F64x2::new1(tail), |s, &ck| s * h + ck);
    let dy = (1..AIRY_TERMS)
        .rev()
        .fold(0.0, |s, k| s * hh + k as f64 * f[k]);
    (y, dy)
}

/// Returns `(e_hi, e_lo)` such that `jn(n, x) ~= exp(e_hi + e_lo)`
/// (`second = false`) or `yn(n, x) ~= -exp(e_hi + e_lo)` (`second = true`),
/// for `n >= ORDER_DEBYE_MIN` and `DEBYE_BELOW_MIN * n <= x < n` with
/// `debye_valid(n, x)`, when `jn_underflows` (`second = false`) or
/// `yn_overflows` (`second = true`) is false.
///
/// With `R = sqrt(n^2 - x^2)` and `w = R / n` (`x = n * sech(α)` and
/// `w = tanh(α)`), Debye's expansions are
///
/// * `jn(n, x) ~= exp(-n * (atanh(w) - w)) / sqrt(2 * π * R) * exp(P(a, b))`
/// * `yn(n, x) ~= -exp(n * (atanh(w) - w)) / sqrt(π * R / 2) * exp(P(-a, -b))`
///
/// where `a = n^2 / R^3` and `b = 1 / R` (see `debye_poly`).
///
/// The absolute error of the exponent is about 2^-67:
/// * `n * (atanh(w) - w)`: `atanh(w)` has a relative error of 2^-100 (see
///   `atanh_precise`), so its error is about `2^-100 * n * w`. With
///   `w <= 2^-6`, it is less than 2^-75 (`n <= 2^31`). With larger `w`,
///   `n * (atanh(w) - w)` is less than 746 (otherwise `jn_underflows` or
///   `yn_overflows` would be true), so
///   `n * w < 746 * w / (atanh(w) - w) < 2^23.2` and the error is less than
///   2^-76. With `w > 0.75`, the relative error of `atanh(w)` is 2^-90, but
///   `n * w < 746 * 0.75 / (atanh(0.75) - 0.75) < 2^11.3`, so the error is
///   less than 2^-78.
/// * Logarithm of `R`: 2^-75 relative to it (see `ln_accurate_parts`), and
///   `ln(R) < 22`.
/// * `P`: the truncation error (less than 2^-73, see `DEBYE_DIST_MIN`), plus
///   the error of the terms of degree 2 and higher, which are less than
///   2^-17.6 and are evaluated with plain `f64` arithmetic (~2^-69).
///
/// The exponential adds 2^-68 (see `exp::Reduced::eval`), so the relative
/// error of the result is about 2^-66.
#[inline(never)]
fn debye_below(n: u32, x: f64, second: bool) -> F64x2 {
    let nf = f64::from(n);
    let r = debye_r(nf, x);
    // w = R / n
    let w = F64x2::div11(r.hi(), nf) + r.lo() / nf;
    // n * (atanh(w) - w)
    let e = (atanh_precise(w) - w) * nf;

    // b = 1 / R and a = n^2 / R^3
    let b = r.recip();
    let y = b * nf;
    let a = b * (y * y);
    let (a, b, e, c) = if second {
        (-a, -b, e, LN_FRAC_PI_2)
    } else {
        (a, b, -e, LN_2PI)
    };
    let p = (DEBYE_A1 * a + DEBYE_B1 * b) + debye_poly(a.hi(), b.hi());

    // ln(2 * π * R) / 2 or ln(π * R / 2) / 2
    let (l_hi, l_lo) = ln_accurate_sum_parts(r.hi(), r.lo());
    let l = (F64x2::add11(l_hi, l_lo) + c) * 0.5;
    (e - l) + p
}

/// Returns `(hi, lo)` such that `hi + lo ~= M(x) * sin(θ(x) + quarter * π/128)`,
/// like `order_large_eval`, where `M` and `θ` are the modulus and the phase
/// of the functions of order `n` (see `debye_above_parts`).
///
/// The relative error is about 2^-65, the one of the sine (see
/// `large_eval`), and the error of the results close to a zero is about
/// 2^-68 times the modulus (the error of the phase).
#[inline(never)]
fn debye_above(n: u32, x: f64, quarter: u8) -> (f64, f64) {
    let (m_hi, m_lo, k, c) = debye_above_parts(n, x);
    // The low part of the sine is not small enough for the product
    let t = c.sin_k(k.wrapping_add(quarter)).normalize();
    mul22(m_hi, m_lo, t.hi(), t.lo())
}

/// Returns `(m_hi, m_lo, k, c)` such that `m_hi + m_lo ~= M(x)` and
/// `θ(x) ~= k * π/128 + c` (modulo `2 * π`), where `M` and `θ` are the
/// modulus and the phase of the functions of order `n`, for
/// `n < x < n^2 * ORDER_LARGE_FACTOR`, with `n >= ORDER_DEBYE_MIN` and
/// `debye_valid(n, x)` or with `debye_far(n, x)`.
///
/// With `R = sqrt(x^2 - n^2)` and `v = R / n` (`x = n * sec(β)` and
/// `v = tan(β)`), Debye's expansions are
///
/// `jn(n, x) + i * yn(n, x) ~= sqrt(2 / (π * R)) * exp(i * θ0) * exp(P(i * a, -i * b))`
///
/// where `θ0 = n * (v - atan(v)) - π/4`, `a = n^2 / R^3` and `b = 1 / R` (see
/// `debye_poly`). So `M(x) = sqrt(2 / (π * R)) * exp(Re(P))` and
/// `θ(x) = θ0 + Im(P)`, where the terms of `P` of even degree are real and the
/// ones of odd degree are imaginary.
///
/// `θ0` is calculated from `n * (v - atan(v))` when `v <= 1`, which is less
/// than `0.22 * n`. Otherwise, `atan(v) = π/2 - atan(1 / v)` and
/// `n * v = R = x - n^2 / (x + R)`, so it is
///
/// `θ0 = x - (2 * n + 1) * π/4 + n * (atan(1 / v) - n / (x + R))`
///
/// where the last term is less than `0.38 * n`. Both are calculated with an
/// error of about `2^-100 * n` (from `atan_precise`), which is less than
/// 2^-69, and they are reduced like `x` (see `reduce_rad`).
///
/// The absolute error of the phase is about 2^-68 (the terms of `P` add the
/// errors of `debye_below`). The relative error of the modulus is about
/// 2^-69, mostly from the evaluation of `Re(P)`, which is less than 2^-17.6,
/// with plain `f64` arithmetic.
#[inline(always)]
fn debye_above_parts(n: u32, x: f64) -> (f64, f64, u8, Reduced) {
    let nf = f64::from(n);
    let r = debye_r(nf, x);
    // b = 1 / R, y = n / R = 1 / v and a = n^2 / R^3
    let b = r.recip();
    let y = b * nf;
    let a = b * (y * y);
    let (p_re, p_im) = debye_poly_complex(a.hi(), b.hi());
    // The terms of degree 1 of P(i * a, -i * b) are imaginary
    let p_im = (DEBYE_A1 * a - DEBYE_B1 * b) + p_im;

    // θ = k * π/128 + c
    let (k, c) = if r.hi() <= nf {
        // θ0 = n * (v - atan(v)) - π/4
        let v = F64x2::div11(r.hi(), nf) + r.lo() / nf;
        let g = (v - atan_precise(v)) * nf;
        let (k, c) = reduce_rad(g.hi());
        let (k, c) = add_phase(k, c, p_im.hi(), p_im.lo() + g.lo());
        (k.wrapping_sub(32), c)
    } else {
        // θ0 = x - (2 * n + 1) * π/4 + n * (atan(y) - n / (x + R))
        let q = F64x2::new1(nf) / (r + x);
        let g = (atan_precise(y) - q) * nf;
        let (k, c) = reduce_rad(g.hi());
        let (k, c) = add_phase(k, c, p_im.hi(), p_im.lo() + g.lo());
        // (2 * n + 1) * π/4 = (64 * n + 32) * π/128 (modulo 2 * π)
        let (k_x, c_x) = reduce_rad(x);
        let k_x = k_x.wrapping_sub((n as u8).wrapping_mul(64).wrapping_add(32));
        // `c_x` and `c` are at most about π/256 (see `add_phase`)
        add_phase(k_x.wrapping_add(k), c_x, c.bh, c.bl)
    };

    // sqrt(2 / (π * R)) * exp(Re(P)), where |Re(P)| < 2^-17.6, so
    // exp(Re(P)) - 1 is approximated by its terms up to the cubic one
    let (w_hi, w_lo) = mul_frac_2_pi(b.hi(), b.lo());
    let w = F64x2::fast_add11(w_hi, w_lo);
    let (m_hi, m_lo) = sqrt_parts(w.hi(), w.lo());
    let e = p_re * (1.0 + p_re * (0.5 + p_re * (1.0 / 6.0)));
    let m = F64x2::fast_add11(m_hi, m_lo + m_hi * e);
    (m.hi(), m.lo(), k, c)
}

// GENERATE: bessel_olver_poly A 1 10 -0.22 0.31
// OLVER_A1[i] is the coefficient of ζ^i
const OLVER_A1: [f64; 11] = [
    f64::from_bits(0xBF723456789ABC3D), // -4.444444444444304e-3
    f64::from_bits(0xBF57FB3BEBF20D9E), // -1.463707463501197e-3
    f64::from_bits(0x3F4725DB8A79BE65), // 7.064172722958782e-4
    f64::from_bits(0x3F460C9721514739), // 6.728876056831596e-4
    f64::from_bits(0x3F242F79DA018B0F), // 1.540027851688265e-4
    f64::from_bits(0xBF0E3B65214E89F5), // -5.7662986513013556e-5
    f64::from_bits(0xBF0A27C8436535E1), // -4.988746552446098e-5
    f64::from_bits(0xBEE5DF8099B197D1), // -1.0429885881357074e-5
    f64::from_bits(0x3ED057A61D4C3931), // 3.8963266504445595e-6
    f64::from_bits(0x3ECA4ABFBB9D9BEF), // 3.13424922364661e-6
    f64::from_bits(0x3E9E29BBCCCC55A5), // 4.4946405068751523e-7
];

// GENERATE: bessel_olver_poly A 2 5 -0.22 0.31
// OLVER_A2[i] is the coefficient of ζ^i
const OLVER_A2: [f64; 6] = [
    f64::from_bits(0x3F46BB7A6E30A55A), // 6.93735867007532e-4
    f64::from_bits(0x3F38294F15688A6F), // 3.686731410274318e-4
    f64::from_bits(0xBF31B1068A55C344), // -2.699509589865654e-4
    f64::from_bits(0xBF37176297388140), // -3.523460102048577e-4
    f64::from_bits(0xBF1A4C5228E6CDB2), // -1.0031939948511903e-4
    f64::from_bits(0x3F12181F097359E6), // 6.902398523620616e-5
];

// GENERATE: bessel_olver_poly B 0 11 -0.22 0.31
// OLVER_B0[i] is the coefficient of ζ^i
const OLVER_B0: [f64; 12] = [
    f64::from_bits(0x3F926E4BDCD91845), // 1.7998872141355333e-2
    f64::from_bits(0x3F823456789ABC9A), // 8.88888888888877e-3
    f64::from_bits(0x3F5AA2A04C9FB099), // 1.62568716268211e-3
    f64::from_bits(0xBF37DFAF9515DDBB), // -3.642848651774462e-4
    f64::from_bits(0xBF33CBBBBFB763B3), // -3.020604487722858e-4
    f64::from_bits(0xBF0EA42A02F32618), // -5.8443576898503126e-5
    f64::from_bits(0x3EF195A00453A3F5), // 1.67698601857592e-5
    f64::from_bits(0x3EEB4C3AC31580D6), // 1.3016591852988231e-5
    f64::from_bits(0x3EC486C4ECEE287C), // 2.4469426089290298e-6
    f64::from_bits(0xBEAA0D2C8005947C), // -7.763940175596619e-7
    f64::from_bits(0xBEA35E2AC74453D6), // -5.772066189758082e-7
    f64::from_bits(0xBE74D7D88CD86267), // -7.764677719573931e-8
];

// GENERATE: bessel_olver_poly B 1 7 -0.22 0.31
// OLVER_B1[i] is the coefficient of ζ^i
const OLVER_B1: [f64; 8] = [
    f64::from_bits(0xBF58756181305A5F), // -1.4928295320310243e-3
    f64::from_bits(0xBF56D71FD15E7994), // -1.394063064686789e-3
    f64::from_bits(0xBF390A7F959A3B16), // -3.820954785646464e-4
    f64::from_bits(0x3F2629B7000A7C2D), // 1.6908987892791505e-4
    f64::from_bits(0x3F266990CF022463), // 1.7099278148105792e-4
    f64::from_bits(0x3F0591C72CC72FEC), // 4.114045222855509e-5
    f64::from_bits(0xBEF231EFB5EF3599), // -1.7352165169450676e-5
    f64::from_bits(0xBEF124FADA2EE4DD), // -1.635022421073584e-5
];

// GENERATE: bessel_olver_poly B 2 3 -0.22 0.31
// OLVER_B2[i] is the coefficient of ζ^i
const OLVER_B2: [f64; 4] = [
    f64::from_bits(0x3F4219127C0966E6), // 5.52305246013432e-4
    f64::from_bits(0x3F4756850D027758), // 7.122182523919647e-4
    f64::from_bits(0x3F2FB5C28589FE96), // 2.4192809611988238e-4
    f64::from_bits(0xBF28C7D68469C449), // -1.8906110399855212e-4
];

/// Returns `jn(n, x)` (`second = false`) or `yn(n, x)` as a sum of two `f64`,
/// for `n >= ORDER_DEBYE_MIN` and `x >= DEBYE_BELOW_MIN * n` where
/// `debye_valid(n, x)` is false, close to `x = n`.
///
/// It uses Olver's uniform asymptotic expansions (DLMF 10.20.4), which are
/// also valid there:
///
/// * `jn(n, x) ~= P * (Ai(t) * (1 + A1 / n^2 + A2 / n^4) + Ai'(t) * (B0 + B1 / n^2 + B2 / n^4) / n^(4/3))`
/// * `yn(n, x) ~= -P * (Bi(t) * (1 + A1 / n^2 + A2 / n^4) + Bi'(t) * (B0 + B1 / n^2 + B2 / n^4) / n^(4/3))`
///
/// where `t = n^(2/3) * ζ`, `P = (4 * ζ / (1 - (x / n)^2))^(1/4) / n^(1/3)`,
/// and `A[k]` and `B[k]` are functions of `ζ` (approximated with the
/// polynomials `OLVER_A1`, ..., `OLVER_B2` for `-0.22 <= ζ <= 0.31`). With
/// `R = sqrt(|n^2 - x^2|)` and `y = R / n`, `ζ` is the solution of
///
/// * `x < n`: `2/3 * ζ^(3/2) = atanh(y) - y`
/// * `x > n`: `2/3 * (-ζ)^(3/2) = y - atan(y)`
///
/// So, with `g = (atanh(y) - y) / y^3` or `(y - atan(y)) / y^3` (see
/// `odd_excess`) and `q = cbrt(3/2 * n * g)`, `t = ±(q * y)^2` and
/// `P = sqrt(2 * q / n)`, which are calculated without cancellations. Their
/// relative errors are about 2^-80.
///
/// Here, `-25.4 < t < 30.9` and `-0.218 < ζ < 0.304` (both extremes with
/// `n = ORDER_DEBYE_MIN`, see `test_olver_args`). The relative error before
/// the final rounding is about 2^-66, relative to the result (`x <= n`) or to
/// the modulus of the functions of order `n` (`x > n`):
/// * Truncation of the expansions: less than 2^-71 (2^-83 with
///   `n >= 4096`).
/// * Approximation of the coefficients: up to 2^-68.5 (`B0` with
///   `n = ORDER_DEBYE_MIN` and the largest `t`).
/// * The error of `t`, about 2^-75 at most, which `Ai(t)` and `Bi(t)`
///   amplify by up to `sqrt(|t|) < 2^2.5`.
/// * `Ai(t)` and `Bi(t)`: about 2^-77 (see `airy`).
/// * The term of the derivatives, which is less than 2^-17 relative to the
///   other one and is calculated with plain `f64` arithmetic: ~2^-68.
#[inline(never)]
fn olver(n: u32, x: f64, second: bool) -> F64x2 {
    let nf = f64::from(n);
    let below = x < nf;
    // y = R / n (`debye_r` requires `x != n`)
    let y = if x == nf {
        F64x2::new1(0.0)
    } else {
        let r = debye_r(nf, x);
        F64x2::div11(r.hi(), nf) + r.lo() / nf
    };
    let q = cbrt_parts(odd_excess(y, below) * (1.5 * nf));
    let qy = q * y;
    let t = if below { qy * qy } else { -(qy * qy) };
    // P = sqrt(2 * q / n)
    let p2 = F64x2::div11(2.0 * q.hi(), nf) + (2.0 * q.lo()) / nf;
    let (p_hi, p_lo) = sqrt_parts(p2.hi(), p2.lo());

    // The terms of the coefficients, with plain `f64` arithmetic
    let n_cbrt = crate::generic::cbrt(nf);
    let zeta = t.hi() / (n_cbrt * n_cbrt);
    let poly = |k: &[f64]| k.iter().rev().fold(0.0, |s, &c| s * zeta + c);
    let e = 1.0 / (nf * nf);
    let sa = (poly(&OLVER_A1) + poly(&OLVER_A2) * e) * e;
    let sb = (poly(&OLVER_B0) + (poly(&OLVER_B1) + poly(&OLVER_B2) * e) * e) / (nf * n_cbrt);

    let (a, da) = airy(t, second);
    let r = F64x2::fast_add11(p_hi, p_lo) * (a + (a.hi() * sa + da * sb));
    if second { -r } else { r }
}

/// Returns `jn(n, x)` rounded to `f64`, for `n >= ORDER_DEBYE_MIN` and finite
/// `DEBYE_BELOW_MIN * n <= x < n^2 * ORDER_LARGE_FACTOR`, when
/// `jn_underflows(n, x, 1075)` is false with `x < n`.
///
/// It uses the Debye expansions (see `debye_below` and `debye_above`), except
/// close to `x = n`, where they are not valid and it uses Olver's (see
/// `olver`).
fn jn_debye(n: u32, x: f64) -> f64 {
    let nf = f64::from(n);
    if !debye_valid(nf, x) {
        olver(n, x, false).to_f64()
    } else if x < nf {
        let e = debye_below(n, x, false);
        let red = ExpReduced::exp_sum(e.hi(), e.lo());
        red.eval_to_f64(exp_poly(red.r))
    } else {
        let (hi, lo) = debye_above(n, x, 64);
        hi + lo
    }
}

/// Returns `yn(n, x)` rounded to `f64`, for `n >= ORDER_DEBYE_MIN` and finite
/// `DEBYE_BELOW_MIN * n <= x < n^2 * ORDER_LARGE_FACTOR`, when
/// `yn_overflows(n, x, 1024)` is false.
///
/// It uses the Debye expansions (see `debye_below` and `debye_above`), except
/// close to `x = n`, where they are not valid and it uses Olver's (see
/// `olver`).
fn yn_debye(n: u32, x: f64) -> f64 {
    let nf = f64::from(n);
    if !debye_valid(nf, x) {
        olver(n, x, true).to_f64()
    } else if x < nf {
        let e = debye_below(n, x, true);
        let red = ExpReduced::exp_sum(e.hi(), e.lo());
        -red.eval_to_f64(exp_poly(red.r))
    } else {
        let (hi, lo) = debye_above(n, x, 0);
        hi + lo
    }
}

/// Returns the sum `hi + lo` as a double-word.
#[inline]
fn normalized((hi, lo): (f64, f64)) -> F64x2 {
    F64x2::add11(hi, lo)
}

/// Returns `true` only if `jn(n, x) < 2^-e`, for `n >= 1` and `0 < x <= n`.
///
/// The backward recurrence of `jn` and the forward one of `yn` take a number
/// of steps proportional to the order also when their result ends up being
/// zero or infinite, so this is checked before them.
///
/// It compares the bound of `jn_log2_bound`, which exceeds the function by a
/// factor of less than about `sqrt(2 * π * n)`, so this returns `true` for
/// all the normal arguments whose result is less than `2^-(e + 18)`.
#[inline]
pub(crate) fn jn_underflows(n: u32, x: f64, e: i32) -> bool {
    // That bound is larger than `(x / n)^n`, and
    // `x / n >= 2^(exponent(x) - ilog2(n) - 1)`, so it does not have to be
    // calculated for most of the arguments (with small orders, only for tiny
    // ones). The exponent of a subnormal `x` is higher here than the actual
    // one, which can only turn the result into `false`.
    let z_exp = i32::from(x.exponent()) - n.ilog2() as i32 - 1;
    if f64::from(n) * f64::from(z_exp) >= -f64::from(e) {
        return false;
    }

    // The error of the calculated bound is much less than 1
    jn_log2_bound(n, x) < -f64::from(e + 1)
}

/// Returns `log2` of an upper bound of `jn(n, x)`, for `n >= 1` and
/// `0 < x <= n`.
///
/// It is Kapteyn's inequality: with `z = x / n` and `w = sqrt(1 - z^2)`,
///
/// `jn(n, x) <= (z * exp(w) / (1 + w))^n`
///
/// The terms of its logarithm nearly cancel when `z` is close to 1, so the
/// absolute error of the result is up to about `n * 2^-49`, plus 2^-51
/// relative to it.
#[inline(never)]
fn jn_log2_bound(n: u32, x: f64) -> f64 {
    let nf = f64::from(n);
    let z = x / nf;
    let w = crate::generic::sqrt((1.0 - z) * (1.0 + z));
    nf * (crate::generic::log2(z / (1.0 + w)) + w * core::f64::consts::LOG2_E)
}

/// Returns `true` only if `yn(n, x) < -2^e`, for `n >= 2` and `x > 0`.
///
/// With `0 < x <= n - 1`, the functions of orders `n - 1` and `n` are
/// positive (first kind) and negative (second kind), so
///
/// `jn(n, x) * yn(n - 1, x) - jn(n - 1, x) * yn(n, x) = 2 / (π * x)`
///
/// implies that `-yn(n, x) > 2 / (π * x * jn(n - 1, x))`, where
/// `π * x / 2 < 2^(exponent(x) + 2)`. So this is `jn_underflows` for the
/// previous order, which returns `true` for all the arguments whose result
/// is less than about `-2^(e + 28)` (the term that is left out is the largest
/// when `x` is close to `n`).
#[inline]
pub(crate) fn yn_overflows(n: u32, x: f64, e: i32) -> bool {
    x <= f64::from(n - 1) && jn_underflows(n - 1, x, e + i32::from(x.exponent()) + 2)
}

/// Returns `jn(n, x)` rounded to `f64`, for `n >= 2` and finite `x > 0`, or
/// zero when `jn_underflows(n, x, e)` is true (so the result is less than
/// 2^-e), for `e <= 1075`. So the methods are only used when
/// `jn_underflows(n, x, 1075)` is false, as they require, and `f32` can use
/// its own limit with a single check.
///
/// For `x < n`, where the function is positive and increasing, the relative
/// error before the final rounding is about 2^-66 (see `jn_backward` and
/// `jn_debye`). Otherwise, where it oscillates, the absolute error is less
/// than about 2^-65 times its modulus (see `jn_forward`, `jn_debye`,
/// `debye_above` and `order_large_eval`), so the results that are small
/// relative to it, close to a zero, are less accurate.
pub(crate) fn jn_finite(n: u32, x: f64, e: i32) -> f64 {
    let nf = f64::from(n);
    if x >= nf * nf * ORDER_LARGE_FACTOR {
        let (hi, lo) = order_large_eval(n, x, 64);
        return hi + lo;
    }

    if x < nf && jn_underflows(n, x, e) {
        return 0.0;
    }
    if n >= ORDER_DEBYE_MIN && x >= DEBYE_BELOW_MIN * nf {
        jn_debye(n, x)
    } else if x < nf {
        jn_backward(n, x)
    } else if debye_far(n, x) {
        let (hi, lo) = debye_above(n, x, 64);
        hi + lo
    } else {
        jn_forward(n, x)
    }
}

/// Returns `true` if `jn(n, x)` and `yn(n, x)` are calculated by
/// `debye_above` for `n < ORDER_DEBYE_MIN` (see `DEBYE_FAR_R_MIN`), which
/// implies `x > n`.
#[inline]
fn debye_far(n: u32, x: f64) -> bool {
    let nf = f64::from(n);
    // The rounding errors only move the limit slightly
    let r2 = (x - nf) * (x + nf);
    (ORDER_DEBYE_FAR_MIN..ORDER_DEBYE_MIN).contains(&n) && r2 >= DEBYE_FAR_R_MIN * DEBYE_FAR_R_MIN
}

/// Returns `jn(n, x)` rounded to `f64`, for `n >= 2` and
/// `n <= x < n^2 * ORDER_LARGE_FACTOR`.
///
/// It uses the recurrence
///
/// `J[k+1](x) = (2 * k / x) * J[k](x) - J[k-1](x)`
///
/// from `j0(x)` and `j1(x)`, with double-word arithmetic. The result is a
/// combination of the functions of both kinds of order `n`, where the one of
/// the second kind comes from the errors of the initial values (2^-68
/// relative to each, which is at most about 2^-66 relative to the modulus)
/// and the rounding errors of each step (about 2^-104). The recurrence does
/// not amplify them because the modulus of the functions of the second kind
/// only starts growing with the order when it is larger than `x`.
fn jn_forward(n: u32, x: f64) -> f64 {
    let t = F64x2::div11(1.0, x);
    let (mut a, mut b) = if x >= J1_LARGE_MIN {
        // Both are calculated from their modulus and phase (`J0_LARGE_MIN` is
        // smaller), which share some parts
        let (a, b) = large_eval_both(x, 64);
        (normalized(a), normalized(b))
    } else {
        (normalized(j0_parts(x)), normalized(j1_parts(x)))
    };
    for k in 1..n {
        let c = (t * f64::from(2 * k)) * b - a;
        a = b;
        b = c;
    }
    b.to_f64()
}

/// Returns the `k` that makes the truncation error of the backward recurrence
/// of `jn_backward` less than about 2^-80 relative to its result, when it
/// starts at the order `m = n + k`, for `n >= 2` and `0 < x < n`.
///
/// That error is about `c / q[k]^2` (see `jn_convergent_terms`), where `c`
/// is at most 1 unless `x` is close to `n`, where the functions of both kinds
/// have similar magnitudes and it grows with the order up to about
/// `n^(2/3) / 8`. So `q[k]` has to reach 2^40 times about the cube root of
/// `n`.
#[inline]
fn jn_backward_terms(n: u32, x: f64) -> u32 {
    // At least 2^39 * n^(1/3)
    jn_convergent_terms(n, x, f64::exp2i_fast(40 + (n.ilog2() / 3) as i16))
}

/// Returns the smallest `k >= 1` such that `q[k] >= q_min`, for `n >= 2` and
/// `0 < x < n`.
///
/// `jn(n, x) / jn(n - 1, x) = 1 / (2 * n / x - 1 / (2 * (n + 1) / x - ...))`,
/// and the denominator of the k-th convergent of this continued fraction,
/// `q[k]`, grows like the functions of the second kind:
///
/// `q[k+1] = (2 * (n + k + 1) / x) * q[k] - q[k-1]`
#[inline]
pub(crate) fn jn_convergent_terms(n: u32, x: f64, q_min: f64) -> u32 {
    let h = 2.0 / x;
    let mut z = f64::from(n) * h;
    let mut q0 = z;
    z += h;
    let mut q1 = q0 * z - 1.0;
    let mut k: u32 = 1;
    while q1 < q_min {
        k += 1;
        z += h;
        (q0, q1) = (q1, z * q1 - q0);
    }
    k
}

/// Returns `jn(n, x)` rounded to `f64`, for `n >= 2` and `0 < x < n`, when
/// `jn_underflows(n, x, 1075)` is false. That also excludes the arguments that
/// are so small that the values below could overflow (`x / n < 2^-539`).
///
/// The recurrence of `jn_forward` is unstable when the order is larger than
/// `x`, so it is used backwards (Miller's algorithm):
///
/// `a[j-1] = (2 * j / x) * a[j] - a[j+1]`
///
/// from `a[m+1] = 0` and `a[m] = 1`, with a large enough `m` (see
/// `jn_backward_terms`), which makes `a[j]` proportional to `J[j](x)` with a
/// relative error that decreases quickly as `j` goes down. The result is
/// `a[n]` normalized with `j0(x)` or `j1(x)` (the one with the largest
/// magnitude).
///
/// The relative error before the final rounding is about 2^-66 (the one of
/// `j0(x)` or `j1(x)`, plus the rounding errors of each step, about 2^-104).
fn jn_backward(n: u32, x: f64) -> f64 {
    let k = jn_backward_terms(n, x);

    // Each step multiplies the magnitude by up to `2 * j / x < 2^556`, so
    // the values are scaled down before they can overflow.
    const SCALE_EXP: i32 = 512;
    let scale_max = f64::exp2i_fast(256);
    let scale_factor = f64::exp2i_fast(-(SCALE_EXP as i16));

    let t = F64x2::div11(1.0, x);
    // a[j + 1] and a[j], and the exponent of their scaling
    let (mut b, mut a) = (F64x2::new1(0.0), F64x2::new1(1.0));
    let mut scale: i32 = 0;
    // a[n] and the exponent of its scaling
    let (mut an, mut an_scale) = (a, scale);

    let m = n + k;
    let mut jf = f64::from(m);
    for j in (1..=m).rev() {
        if j == n {
            (an, an_scale) = (a, scale);
        }
        let c = (t * (2.0 * jf)) * a - b;
        b = a;
        a = c;
        jf -= 1.0;
        if a.hi().abs() > scale_max {
            a *= scale_factor;
            b *= scale_factor;
            scale += SCALE_EXP;
        }
    }

    // `a` and `b` are proportional to `j0(x)` and `j1(x)`
    let r = if a.hi().abs() >= b.hi().abs() {
        normalized(j0_parts(x)) * an / a
    } else {
        normalized(j1_parts(x)) * an / b
    };
    // The result can be subnormal, so it is rounded only once
    r.scalbn_to_f64(an_scale - scale)
}

/// Returns `yn(n, x)` rounded to `f64`, for `n >= 2` and finite `x > 0`, or
/// negative infinity when `yn_overflows(n, x, e)` is true (so the result is
/// less than -2^e), for `e <= 1024` (like `jn_finite`).
///
/// It uses the recurrence of `jn_forward` from `y0(x)` and `y1(x)`, which is
/// stable for the functions of the second kind with any order. The exceptions
/// are `x >= n^2 * ORDER_LARGE_FACTOR` (see `order_large_eval`), the orders
/// from `ORDER_DEBYE_MIN` with `x >= DEBYE_BELOW_MIN * n`, where the
/// asymptotic expansions of `yn_debye` are faster, and the arguments far
/// beyond smaller orders (see `debye_far`).
///
/// The absolute error before the final rounding is less than about 2^-65
/// times the modulus of the functions of order `n`, like in `jn_forward`
/// (with all the methods). So the relative error is as small for `x < n`,
/// where the function is negative and its magnitude is the modulus, while the
/// results that are small relative to the modulus, close to a zero, are less
/// accurate.
pub(crate) fn yn_finite(n: u32, x: f64, e: i32) -> f64 {
    let nf = f64::from(n);
    if x >= nf * nf * ORDER_LARGE_FACTOR {
        let (hi, lo) = order_large_eval(n, x, 0);
        return hi + lo;
    }

    if yn_overflows(n, x, e) {
        return f64::NEG_INFINITY;
    }
    if n >= ORDER_DEBYE_MIN && x >= DEBYE_BELOW_MIN * nf {
        return yn_debye(n, x);
    }
    if debye_far(n, x) {
        let (hi, lo) = debye_above(n, x, 0);
        return hi + lo;
    }

    // It can overflow with subnormal `x`, and so does the result
    let t = F64x2::div11(1.0, x);
    let (mut a, mut b) = if x >= Y0_LARGE_MIN {
        // Both are calculated from their modulus and phase (`Y1_LARGE_MIN` is
        // smaller), which share some parts
        let (a, b) = large_eval_both(x, 0);
        (normalized(a), normalized(b))
    } else {
        (normalized(y0_parts(x)), normalized(y1_parts(x)))
    };

    // The product of each step is larger in magnitude than its result, so
    // it could overflow when the result does not (with large orders, that
    // would be the results up to several percent below the largest finite
    // value). So the values are scaled down once their magnitude gets large,
    // which only happens with `x` smaller than the order, where it keeps
    // growing with it. Before that, a magnitude of at most 2^256 implies
    // that `1 / x` is less than 2^257, so the products are less than
    // 2^32 * 2^257 * 2^256.
    const SCALE_EXP: i32 = 512;
    let scale_min = f64::exp2i_fast(256);
    let scale_factor = f64::exp2i_fast(-(SCALE_EXP as i16));
    let mut scale: i32 = 0;
    for k in 1..n {
        if !b.hi().is_finite() {
            break;
        }
        if scale == 0 && b.hi().abs() > scale_min {
            a *= scale_factor;
            b *= scale_factor;
            scale = SCALE_EXP;
        }
        let c = (t * f64::from(2 * k)) * b - a;
        a = b;
        b = c;
    }
    if b.hi().is_finite() {
        // It overflows here when the result does
        b.scalbn_to_f64(scale)
    } else {
        // The scaled values (or `1 / x`) only overflow when the result is
        // far beyond the largest finite value
        f64::NEG_INFINITY
    }
}

/// Smallest argument of `AIRY_TBL` (see `airy`).
const AIRY_T_MIN: f64 = -26.0;

// The arguments of `airy` in `olver` are at most 1/8 beyond the first and
// the last points of `AIRY_TBL`
const _: () = {
    let t_max = AIRY_T_MIN + (AIRY_TBL.len() - 1) as f64 / 4.0;
    assert!(AIRY_T_MIN - 0.125 < -25.4 && 30.9 < t_max + 0.125);
};

// GENERATE: bessel_airy_table -26 31 4
// AIRY_TBL[i] = [Ai(t), Ai'(t), Bi(t), Bi'(t)], with t = -26 + i / 4
static AIRY_TBL: [[F64x2; 4]; 229] = [
    // t = -26
    [
        F64x2::from_bits(0x3FCDD79FD9C6EF65, 0xBC40E8198FD50075), // 2.331428349842738205188307591097e-1
        F64x2::from_bits(0xBFDD2C03359B8B14, 0xBC77ECBDD1DE6E40), // -4.558113120566969141869458593856e-1
        F64x2::from_bits(0x3FB6FF2654BC7D51, 0xBC52925B5FE280DC), // 8.983077591664324498567486816856e-2
        F64x2::from_bits(0x3FF308E7C100B308, 0xBC9374401AEB06DE), // 1.189674142758294127737359161873e0
    ],
    // t = -25.75
    [
        F64x2::from_bits(0xBF9199420E4C9603, 0xBC387F92A74A9FFF), // -1.718619549060030956689492431332e-2
        F64x2::from_bits(0xBFF44A238A6B7CD8, 0xBC930E854192ED34), // -1.268100300505286947619079757966e0
        F64x2::from_bits(0x3FCFFB8C0721A231, 0x3C6977F77B7344E3), // 2.498641047307202211465778628144e-1
        F64x2::from_bits(0xBFB5B4805D9456DF, 0xBC5579D54AB7039C), // -8.478548321394495845217268780267e-2
    ],
    // t = -25.5
    [
        F64x2::from_bits(0xBFCF3DC434AF04D3, 0xBC6CAF6219B06053), // -2.440724618191213293240266282497e-1
        F64x2::from_bits(0xBFD32BD653F25BF3, 0x3C6E7C052FCEB524), // -2.995506114761489625442193454037e-1
        F64x2::from_bits(0x3FAE20FE787CD23F, 0xBC27DF0D797B49F6), // 5.884547444148723679103141835674e-2
        F64x2::from_bits(0xBFF3B606F9583000, 0x3C9D5F86824CE066), // -1.231940244676479676809630773313e0
    ],
    // t = -25.25
    [
        F64x2::from_bits(0xBFC0C8D9141D98A2, 0xBC4970C71B6CDA53), // -1.311293934078934825087387495906e-1
        F64x2::from_bits(0x3FF14051BBBFCC4E, 0x3C8F855A3AD8C329), // 1.078202947039943958852846545977e0
        F64x2::from_bits(0xBFCB7F70CE7A8957, 0xBC481972E2197C03), // -2.148266800114224410228146805160e-1
        F64x2::from_bits(0xBFE52752FCFA862E, 0xBC79BEF07DD78E86), // -6.610503140438533044175871326019e-1
    ],
    // t = -25
    [
        F64x2::from_bits(0x3FC4EE705E0D3E32, 0xBC6A6E11C438A058), // 1.635265788304294694863710011142e-1
        F64x2::from_bits(0x3FEECBCEBBA24F80, 0x3C29E66CEFA5E29D), // 9.623788513876974100384205686412e-1
        F64x2::from_bits(0xBFC8984450B5D9EF, 0x3C60ED032EF4B7A8), // -1.921468156903780237647036128162e-1
        F64x2::from_bits(0x3FEA1A603BBBCA54, 0xBC8251D36082A901), // 8.157197157546058578776189286478e-1
    ],
    // t = -24.75
    [
        F64x2::from_bits(0x3FCE0DA0B90D57FB, 0xBC681960F21FDD09), // 2.347908881735632675561156028125e-1
        F64x2::from_bits(0xBFDDCF75200A8DED, 0xBC7693B95CF13A88), // -4.657872021582239648188817884936e-1
        F64x2::from_bits(0x3FB81719659B125F, 0x3C52BF98CEC22423), // 9.410246592148884388273825977324e-2
        F64x2::from_bits(0x3FF2B45B58A6C289, 0xBC987F1B21D909B1), // 1.169032427115015119687872220708e0
    ],
    // t = -24.5
    [
        F64x2::from_bits(0xBF8A78F85A22ED7B, 0x3C26AD88DFD541D3), // -1.292604470324109259698038793893e-2
        F64x2::from_bits(0xBFF40F39FEFFE863, 0xBC980E009CE1AD9C), // -1.253717418758719087857808792660e0
        F64x2::from_bits(0x3FD03568BA17B250, 0xBC69DC4E74F20B5F), // 2.532598321256829153582823561069e-1
        F64x2::from_bits(0xBFAF6F74C1077D21, 0x3C47B685E021F9CA), // -6.139721733392833772833182514897e-2
    ],
    // t = -24.25
    [
        F64x2::from_bits(0xBFCF39B1058FE63F, 0x3C54CC398BE81279), // -2.439481045874334019107219635442e-1
        F64x2::from_bits(0xBFD6BAA5D4F0361D, 0x3C6A6E95B7D62757), // -3.551420764423214745109224419777e-1
        F64x2::from_bits(0x3FB254D5A995F622, 0xBC427B13979F06AD), // 7.160697355535508367734013210341e-2
        F64x2::from_bits(0xBFF33593776CC2DE, 0xBC9679A430422E93), // -1.200580088141741075519064340654e0
    ],
    // t = -24
    [
        F64x2::from_bits(0xBFC32DD86B12AB90, 0x3C6A98553BF84146), // -1.498365900818865332891997815711e-1
        F64x2::from_bits(0x3FF02387C694925C, 0x3C8C6CEA28FEAF79), // 1.008674407677197009574901476409e0
        F64x2::from_bits(0xBFCA651FA1613666, 0xBC6D1C7F4406023E), // -2.062110460268968396797630539635e-1
        F64x2::from_bits(0xBFE78EF8ACCB1E0B, 0x3C7BBD6AADD62A13), // -7.362025618430732868715483219987e-1
    ],
    // t = -23.75
    [
        F64x2::from_bits(0x3FC24970154569D8, 0x3C34E7D46DA4B598), // 1.428661445487040560173636336915e-1
        F64x2::from_bits(0x3FF08C29FFC819A6, 0x3C94755BEB1D5F2C), // 1.034219741007839825566844891204e0
        F64x2::from_bits(0xBFCB1FBF8CAC308B, 0xBC6A143506357051), // -2.119063794019705295602436249124e-1
        F64x2::from_bits(0x3FE6356BD5964B4F, 0xBC8556E08CDB8CA4), // 6.940211459753785246476409228516e-1
    ],
    // t = -23.5
    [
        F64x2::from_bits(0x3FCFE0D5D21D80E1, 0x3C6EF31C7C243D2C), // 2.490489269212110968628830028005e-1
        F64x2::from_bits(0xBFD289FCB2B49161, 0x3C7DE2DB75A6CB44), // -2.896720643205039931749451258071e-1
        F64x2::from_bits(0x3FAEDFB70C39FA9A, 0x3C402E1994055882), // 6.030056022881531871969254190423e-2
        F64x2::from_bits(0x3FF353D393B1AB6A, 0xBC950CC2DCE98696), // 1.207965447358868276394287906178e0
    ],
    // t = -23.25
    [
        F64x2::from_bits(0x3FA0573195F16F12, 0x3C3B5AC157B01A2D), // 3.191523509369011283939250767784e-2
        F64x2::from_bits(0xBFF3A9CB86151083, 0xBC93CE4C66489E4F), // -1.228953861006659333702563609813e0
        F64x2::from_bits(0x3FD050F4B4B68730, 0x3C7663A94EBD6A64), // 2.549411549353975042682599905713e-1
        F64x2::from_bits(0x3FC40C8A5587A629, 0xBC583D9755D9334E), // 1.566327016474116321126493810256e-1
    ],
    // t = -23
    [
        F64x2::from_bits(0xBFCD0C2CD0CB925F, 0xBC4854D3377B3639), // -2.269340533740828751164777471195e-1
        F64x2::from_bits(0xBFE2CB72F57128CF, 0xBC72E100F6EFA8FD), // -5.873350900449397946058793869840e-1
        F64x2::from_bits(0x3FBF383E861C16A0, 0x3C34061862D72CE5), // 1.219519689003285567193980190447e-1
        F64x2::from_bits(0xBFF164755F19F6F5, 0xBC83DDBEC29914E2), // -1.087025996675433682271322310383e0
    ],
    // t = -22.75
    [
        F64x2::from_bits(0xBFC93C30131CC50E, 0xBC22B5027539EF26), // -1.971492856332947579208920105517e-1
        F64x2::from_bits(0x3FE96908B4F2D4A7, 0xBC8BD7B4D60BA30E), // 7.940715345527181240291943705329e-1
        F64x2::from_bits(0xBFC55E1CE2B104D1, 0x3C6E166FD56B92AB), // -1.669345957653036874279873596251e-1
        F64x2::from_bits(0xBFEE2669690C677C, 0xBC81A7790025F1A6), // -9.421889354446047222076643245914e-1
    ],
    // t = -22.5
    [
        F64x2::from_bits(0x3FB4EBD4D3B28CFD, 0x3C5296185DE085C9), // 8.172350090403664572018375846238e-2
        F64x2::from_bits(0x3FF2ABC5240C28E5, 0x3C6B5C74B7653FFF), // 1.166936055002730355454278513853e0
        F64x2::from_bits(0xBFCF76EDE75C8015, 0xBC4FA0EAA47464E8), // -2.458169345563880093682757606844e-1
        F64x2::from_bits(0x3FD8A2936BB276A9, 0x3C6B8DB23317BD4D), // 3.849228431443899153964039761595e-1
    ],
    // t = -22.25
    [
        F64x2::from_bits(0x3FD09673A60B210F, 0x3C5B48013AFBD75D), // 2.591828462141316295234947459020e-1
        F64x2::from_bits(0x3FB5D468EA21EE38, 0x3C5B75E8FC1338F4), // 8.527236667753690975523702059754e-2
        F64x2::from_bits(0xBF91E111A6033BCD, 0xBC204B975D881752), // -1.746013236771286727714740978685e-2
        F64x2::from_bits(0x3FF38EE2BF80597C, 0x3C98F0349F4E7EFE), // 1.222384212543686652288284751215e0
    ],
    // t = -22
    [
        F64x2::from_bits(0x3FBDBB9F8E05E390, 0x3C58E970657615B3), // 1.161441537605141292761599777811e-1
        F64x2::from_bits(0xBFF17A85CB0710D9, 0x3C81FCE2C8650B42), // -1.092412751270833672794243906391e0
        F64x2::from_bits(0x3FCDD8E1072BC8D3, 0xBC60593FF6FF78D3), // 2.331811223044214160957826789264e-1
        F64x2::from_bits(0x3FE1847AF59F8266, 0xBC8D845B6F525E91), // 5.474219128919456987865732164697e-1
    ],
    // t = -21.75
    [
        F64x2::from_bits(0xBFC5B9B876329F46, 0xBC69582535142687), // -1.697302415025154330720226400717e-1
        F64x2::from_bits(0xBFEDB3BCC7A27EAA, 0xBC8BBC0B5B85DD03), // -9.281906031719320752495857825522e-1
        F64x2::from_bits(0x3FC96BD83AFD850D, 0xBC5DD466B1FD62F1), // 1.986036575354223278887392025452e-1
        F64x2::from_bits(0xBFE941EF1503B697, 0xBC8570F4F3420E37), // -7.892985735358761062462715979432e-1
    ],
    // t = -21.5
    [
        F64x2::from_bits(0xBFD0064EFCC27C75, 0x3C4EA87E9AA3765B), // -2.503850429874952510167906218805e-1
        F64x2::from_bits(0x3FD6B6B91C4F0D2A, 0x3C4AA0B6E4BCC196), // 3.549025322392532795227839557242e-1
        F64x2::from_bits(0xBFB3C1355C839104, 0x3C5867ED607BAAC8), // -7.716687687531381149376152670163e-2
        F64x2::from_bits(0xBFF29727BDA06553, 0x3C8A87EEE481728D), // -1.161903134087272258573675095755e0
    ],
    // t = -21.25
    [
        F64x2::from_bits(0xBF9F246833740A82, 0xBBF22BD2ED5023C8), // -3.041231931730204596836981981067e-2
        F64x2::from_bits(0x3FF33EDDF046B446, 0x3C65B8CB91200289), // 1.202848375872478846032698287510e0
        F64x2::from_bits(0xBFD0B45A189CAEA4, 0x3C60262C3558DD44), // -2.610078087188674108484777159244e-1
        F64x2::from_bits(0xBFC2568EBCB870D6, 0x3C6497CFF7A8CAA6), // -1.432665273176840136713535742052e-1
    ],
    // t = -21
    [
        F64x2::from_bits(0x3FCCF950ABC2F44C, 0xBC51C85A0E2212EA), // 2.263584936789889662543187309113e-1
        F64x2::from_bits(0x3FE3E1A4E5F171E3, 0xBC8F42940FAF3FA8), // 6.212944499089270724080801422854e-1
        F64x2::from_bits(0xBFC14743976ABBEE, 0x3C61026221ABB71C), // -1.349873055782774759544943384950e-1
        F64x2::from_bits(0x3FF0924A849217AC, 0xBC944BE71E6B5808), // 1.035715597013658602702448398977e0
    ],
    // t = -20.75
    [
        F64x2::from_bits(0x3FCBD6D91D0A5B18, 0xBC61F9135195947E), // 2.174941436575969977199784865920e-1
        F64x2::from_bits(0xBFE5D12DC72812FB, 0x3C82736DC64644A8), // -6.817845239425815675132530323115e-1
        F64x2::from_bits(0x3FC33B2FA27571BF, 0xBC5DBADEAA9C402F), // 1.502437155967423989470690015533e-1
        F64x2::from_bits(0x3FEFC30DC1A9506B, 0xBC48E0DA471C095E), // 9.925602705429928924401878043416e-1
    ],
    // t = -20.5
    [
        F64x2::from_bits(0xBFA6D92D5BBBDC36, 0x3C2416616205320E), // -4.462568039701190981636756810425e-2
        F64x2::from_bits(0xBFF2F163C0049C30, 0xBC97652E422DDA5B), // -1.183933019705147497002823901988e0
        F64x2::from_bits(0x3FD0BA2842A3E737, 0xBC34957A25AA008E), // 2.613621378692302921142054369052e-1
        F64x2::from_bits(0xBFC97481EE495D13, 0x3C614625BD3520B4), // -1.988680280216859810908253503544e-1
    ],
    // t = -20.25
    [
        F64x2::from_bits(0xBFD063DDC48EE2AF, 0x3C5B07E7D9E6B8FF), // -2.560953540337704223090253313493e-1
        F64x2::from_bits(0xBFD4DE8CC91EDE46, 0xBC72BEAC73C517D0), // -3.260833705967097924598327719615e-1
        F64x2::from_bits(0x3FB25ECE0AB5B1A6, 0x3C55A8A528CD6FEB), // 7.175910722225484675876721213138e-2
        F64x2::from_bits(0xBFF26CCF5F95BA9C, 0x3C6BC8F728DE8921), // -1.151564954164917319475143245258e0
    ],
    // t = -20
    [
        F64x2::from_bits(0xBFC69479D94E9662, 0x3C633F42AF1C6953), // -1.764061270779846895901922922195e-1
        F64x2::from_bits(0x3FEC9255202FE181, 0xBC73946820C67673), // 8.928628567364712383984099341143e-1
        F64x2::from_bits(0xBFC99E2A3617A808, 0xBC6B0CBAC8A80452), // -2.001393093226513492836047993359e-1
        F64x2::from_bits(0xBFE95362FB2E4310, 0xBC842C57EE4CBC20), // -7.914290338395364793562815314625e-1
    ],
    // t = -19.75
    [
        F64x2::from_bits(0x3FBA2E465D081917, 0x3C47571C765B6543), // 1.022685982826242343294565233411e-1
        F64x2::from_bits(0x3FF19B48C55649BE, 0xBC97100ED79A48E8), // 1.100411196563626434673488240058e0
        F64x2::from_bits(0xBFCFA8098DE99AAC, 0xBC54C8B5A9A570A7), // -2.473155920917710127275627036730e-1
        F64x2::from_bits(0x3FDCE3418DC5BC1C, 0x3C7A089C9E236413), // 4.513706097122296649213136668953e-1
    ],
    // t = -19.5
    [
        F64x2::from_bits(0x3FD123A3C0A2978A, 0xBC3AD0A754CFA9BB), // 2.678002721025839457554509636200e-1
        F64x2::from_bits(0x3FB6763330F3AEA7, 0xBC5FDB488A92E6B4), // 8.774108834375713570132023758930e-2
        F64x2::from_bits(0xBF938CC26C7E2E3A, 0x3C3BA0197D0E8EDD), // -1.909164225757662414329380469194e-2
        F64x2::from_bits(0x3FF2EAEC3108E970, 0xBC90AAABEB35594D), // 1.182354156057581149226861595555e0
    ],
    // t = -19.25
    [
        F64x2::from_bits(0x3FC1C53DEA3B3112, 0x3C6B103BA7393D82), // 1.388318437296978797542086596438e-1
        F64x2::from_bits(0xBFF02CA3A6EFD85F, 0xBC9CBB1BFFF496FC), // -1.010898258301970875399133897300e0
        F64x2::from_bits(0x3FCD8B3671B34721, 0x3C44813FD21CD79D), // 2.308109336068407146962729815196e-1
        F64x2::from_bits(0x3FE396976098A5F4, 0x3C77F97BCE99EAFB), // 6.121327292104852476951766140289e-1
    ],
    // t = -19
    [
        F64x2::from_bits(0xBFC221F4EBA7CF61, 0xBC4500E94E39FACD), // -1.416612768804226563669522747023e-1
        F64x2::from_bits(0xBFF014521DDA6005, 0xBC7F951345B6D32B), // -1.004961125005139593506444282180e0
        F64x2::from_bits(0x3FCD749B9E09BE82, 0x3C5694A4EF87A63D), // 2.301210900945883146676358262322e-1
        F64x2::from_bits(0xBFE3A9C4DCAFBBEA, 0x3C623EF5DB47EBF4), // -6.144737539560740567619935164991e-1
    ],
    // t = -18.75
    [
        F64x2::from_bits(0xBFD150187D1464B9, 0x3C75A4314DC4727C), // -2.705136510354111768402440539935e-1
        F64x2::from_bits(0x3FB33FCC72C1D345, 0x3C3AAB9BAE82216F), // 7.519223977744811527020074971174e-2
        F64x2::from_bits(0xBF92A25181F35D34, 0x3C332532519F7019), // -1.819732051521079259656204452298e-2
        F64x2::from_bits(0xBFF2BEFDDBB4B22C, 0x3C96EF7304D5AA17), // -1.171628816820098704462077558354e0
    ],
    // t = -18.5
    [
        F64x2::from_bits(0xBFBCB1D5A497A3B2, 0x3C51CFF12C8E0A30), // -1.120885397755404761189093926743e-1
        F64x2::from_bits(0x3FF108C81B7E00D7, 0x3C962001A40B5E1D), // 1.064643962279708443597562571989e0
        F64x2::from_bits(0xBFCFBA3A0F23FE09, 0xBC6D87080B4A5446), // -2.478706907204329252072023669723e-1
        F64x2::from_bits(0xBFDF11F94FF9973D, 0xBC7EC4249CD88676), // -4.854720383649356585787704036974e-1
    ],
    // t = -18.25
    [
        F64x2::from_bits(0x3FC50F9A5E9DB4CB, 0x3C67E70712DAA204), // 1.645386659614927852689438720621e-1
        F64x2::from_bits(0x3FEDD8C1A9039AA5, 0x3C85C28378BF7B05), // 9.327095318847865956787452470149e-1
        F64x2::from_bits(0xBFCBE0CB362D3C2E, 0x3C5C1E99B0C4706B), // -2.177976622137704426717161687932e-1
        F64x2::from_bits(0x3FE665F1F61D5E0A, 0xBC76114A8735AEDD), // 6.999444777420091665206887923706e-1
    ],
    // t = -18
    [
        F64x2::from_bits(0x3FD15B6A4A5206A9, 0x3C7DA713CC9BA2BE), // 2.712045408044142215843503642002e-1
        F64x2::from_bits(0xBFC45B631DCC5F2F, 0x3C41CA6A298EC953), // -1.590389152049680157736995756187e-1
        F64x2::from_bits(0x3FA3A58F0E6B569A, 0x3C121A0C41DD5851), // 3.837248850838399807480063842072e-2
        F64x2::from_bits(0x3FF26B43288C61DB, 0x3C9435E745F05986), // 1.151187094108641798715365191767e0
    ],
    // t = -17.75
    [
        F64x2::from_bits(0x3FB9A8669E98D592, 0x3C5FC88E77CEFD48), // 1.002258431625107428871723528974e-1
        F64x2::from_bits(0xBFF13B0B6C7FADA1, 0x3C6BA89ADF83B979), // -1.076915191481795585288610248705e0
        F64x2::from_bits(0x3FD061543F6D17DC, 0xBC76EB2AC6FD66BF), // 2.559404963870954642311283598438e-1
        F64x2::from_bits(0x3FDB4189125C2B7D, 0x3C23A22A76580C96), // 4.258749656007053135723744772913e-1
    ],
    // t = -17.5
    [
        F64x2::from_bits(0xBFC619BE0319E520, 0x3C63FFACCF6F9FE3), // -1.726605906622262678196076382759e-1
        F64x2::from_bits(0xBFECE08048A6DEC9, 0xBC8AB097DB2456C9), // -9.024049204808416898591519658455e-1
        F64x2::from_bits(0x3FCB890F69BBE6EF, 0x3C59DD266BF6D5D2), // 2.151202455786953349775998676590e-1
        F64x2::from_bits(0xBFE7040291FBBF8F, 0x3C710552867E4E08), // -7.192395068395728138260863192857e-1
    ],
    // t = -17.25
    [
        F64x2::from_bits(0xBFD18650400D9635, 0xBC70ADA293BEF3D5), // -2.738228440779012096471991233406e-1
        F64x2::from_bits(0x3FC524FFD6CB884C, 0xBC24F1C0C72B4F95), // 1.651916312030742289714428756260e-1
        F64x2::from_bits(0xBFA4DA3DB8930F35, 0x3C448FF7E2E0C3FD), // -4.072754743677223703496920161489e-2
        F64x2::from_bits(0xBFF234D2E2899F0A, 0xBC9242D9B5E6F254), // -1.137896428024705538579791975480e0
    ],
    // t = -17
    [
        F64x2::from_bits(0xBFBAF278594059E7, 0x3BE5EDAD2770C2CA), // -1.052623002909523902321130354513e-1
        F64x2::from_bits(0x3FF0F05F3D177EA6, 0x3C9BF7E7EDEA1F77), // 1.058684576644660077394894917170e0
        F64x2::from_bits(0xBFD074EA38D53DE8, 0xBC6A53CF3BB04A74), // -2.571359210023431821447261248492e-1
        F64x2::from_bits(0xBFDC04F2F4AC62D3, 0x3C5CEC0A4F34A973), // -4.378020657909875094650659795462e-1
    ],
    // t = -16.75
    [
        F64x2::from_bits(0x3FC5460626B8E05C, 0xBC5BD2CC75890E49), // 1.661994637559044358418435018701e-1
        F64x2::from_bits(0x3FED68C31555FA0B, 0xBC825D64AE2B9D78), // 9.190383354869607571757368804689e-1
        F64x2::from_bits(0xBFCCAA2D3DFBC786, 0xBC68183C1B220A6C), // -2.239433815626023609590459810341e-1
        F64x2::from_bits(0x3FE5A90347DE4DF3, 0xBC8998D19290595B), // 6.768814471892582007636829675971e-1
    ],
    // t = -16.5
    [
        F64x2::from_bits(0x3FD1D8FB2EF92345, 0xBC764C87BDCBB152), // 2.788684805605508383085624965939e-1
        F64x2::from_bits(0xBFB8392F76681165, 0x3C37CA7FA0364282), // -9.462257996353213999977071365139e-2
        F64x2::from_bits(0x3F98EAE825D985E7, 0xBC21CD7BEF0CB135), // 2.433359843269569271524904183766e-2
        F64x2::from_bits(0x3FF2217E51B888CF, 0x3C5376E3B93549AB), // 1.133177108022710437353500783401e0
    ],
    // t = -16.25
    [
        F64x2::from_bits(0x3FC05B15674E577B, 0xBC62647DAA47F53B), // 1.277796511100709107193885556984e-1
        F64x2::from_bits(0xBFF01C58FF7C116D, 0x3C9D5684C1845401), // -1.006920812594297785301639143939e0
        F64x2::from_bits(0x3FD0045680C3F476, 0x3C7B3B36E36D91BD), // 2.502647645449721750248817459448e-1
        F64x2::from_bits(0x3FE09B5B2A7ADF3A, 0x3C828655750D1992), // 5.189643697440040936374001133837e-1
    ],
    // t = -16
    [
        F64x2::from_bits(0xBFC24FB8E8F9E13A, 0xBC5C747BB970CC15), // -1.430579316690996977774693643023e-1
        F64x2::from_bits(0xBFEF314532C2256B, 0x3C75E889AD5DEC9E), // -9.747644416212727179572594096248e-1
        F64x2::from_bits(0x3FCF1EA8D02470C6, 0x3C4FA7EE920FEBF8), // 2.431231514282272166875883325401e-1
        F64x2::from_bits(0xBFE230C9CF9CA6D2, 0x3C8E3F81A9F69774), // -5.684556059761353727226430431125e-1
    ],
    // t = -15.75
    [
        F64x2::from_bits(0xBFD21BBA90249E39, 0xBC78E0B1288DE39D), // -2.829424293181293730677979191054e-1
        F64x2::from_bits(0xBFAAEE0DDF332CF1, 0xBC2A464F04AE68B6), // -5.259745931251725502727587737747e-2
        F64x2::from_bits(0x3F88D31FC89D0BEE, 0xBC2E24D5B7E9F015), // 1.212143736814216403615027800339e-2
        F64x2::from_bits(0xBFF1F6C3F2DF0091, 0xBC82F57CEDA2C397), // -1.122745465007060771815239175258e0
    ],
    // t = -15.5
    [
        F64x2::from_bits(0xBFC54E2AA3A6F89F, 0xBC4E7A2899E5B626), // -1.664479540904197673881618281096e-1
        F64x2::from_bits(0x3FECF54066B2A403, 0x3C80D6C7B224A904), // 9.049379354302121995067425701019e-1
        F64x2::from_bits(0xBFCD81E4B339D579, 0x3C5B57E7B8E98EA2), // -2.305265307547122107351091693995e-1
        F64x2::from_bits(0xBFE516F2082AD467, 0xBC5448E19D5B5380), // -6.590509566800734119883009851119e-1
    ],
    // t = -15.25
    [
        F64x2::from_bits(0x3FB966A4A35C3AE6, 0xBC5C3E9D91BAA9E0), // 9.922245968139583536635557440744e-2
        F64x2::from_bits(0x3FF0C0C7DD27D3C7, 0xBC90D3F8F258A936), // 1.047065605057683587611412017409e0
        F64x2::from_bits(0xBFD121F60EB24B2C, 0xBC4156317AE0F9BE), // -2.676978248138535151279664175117e-1
        F64x2::from_bits(0x3FD884CE40794CA5, 0xBC63E909061D4351), // 3.831058148982243085873171976093e-1
    ],
    // t = -15
    [
        F64x2::from_bits(0x3FD1CE50BC1DC8EE, 0x3C6453703831F307), // 2.782174908708289295276215087712e-1
        F64x2::from_bits(0x3FD16E9436F1E59E, 0x3C6A25763B13AE78), // 2.723742043086420208257838802340e-1
        F64x2::from_bits(0xBFB1B247CECB6704, 0xBC5CEF8C5ABC0EA3), // -6.912659453101006118592810401063e-2
        F64x2::from_bits(0x3FF1390E679F0519, 0x3C93A461EEB27EB8), // 1.076429753084374786744196080003e0
    ],
    // t = -14.75
    [
        F64x2::from_bits(0x3FCBAF9B1238D108, 0x3C4F103DB3138F75), // 2.162965620970072089618732092633e-1
        F64x2::from_bits(0xBFE73B6DF33C4820, 0x3C83ADFE3BE73943), // -7.260045767491227422524628765598e-1
        F64x2::from_bits(0x3FC851482DDB4F57, 0x3C47261A4E7A79EB), // 1.899805282506530104809403518929e-1
        F64x2::from_bits(0x3FEAAFD175F4CE69, 0xBC872FDFF5F555D0), // 8.339621833642470602690978383537e-1
    ],
    // t = -14.5
    [
        F64x2::from_bits(0xBF9F54EE06BFE1B0, 0x3C35BD8D7B00B823), // -3.059741893955142282119372095576e-2
        F64x2::from_bits(0xBFF1866F995A28FE, 0x3C8CD6950EB047E5), // -1.095321272880539215033628254669e0
        F64x2::from_bits(0x3FD26645DDF0BB98, 0xBC7B777C8FB60D0A), // 2.874922435175277552383898228014e-1
        F64x2::from_bits(0xBFBC8F5782AED737, 0xBC556777C0C290E1), // -1.115622228670333134564775560474e-1
    ],
    // t = -14.25
    [
        F64x2::from_bits(0xBFD0280EE7314F1F, 0x3C5A8B0EEC19A8C9), // -2.524449594250074922105819210614e-1
        F64x2::from_bits(0xBFE179CE9FA5F2BA, 0xBC705393112D6C2D), // -5.461190336351002958037982036648e-1
        F64x2::from_bits(0x3FC25DE04C079D76, 0x3C6BE15BF34C5C16), // 1.434898730505939755572273562914e-1
        F64x2::from_bits(0xBFEE6A71B23AD859, 0x3C76576D66247A4F), // -9.504936677028680087017031093710e-1
    ],
    // t = -14
    [
        F64x2::from_bits(0xBFD105DF95F55CDE, 0xBC4615077C1770BA), // -2.659834827840777983847914327271e-1
        F64x2::from_bits(0x3FDC5A850382AA2A, 0xBC6C15605CE8D5B3), // 4.430248770028436411715168920642e-1
        F64x2::from_bits(0xBFBEA266D3B93B23, 0xBC5EBC5797A4E62E), // -1.196655527976245231326009740984e-1
        F64x2::from_bits(0xBFEFEACC30E10016, 0xBC74819573F34C58), // -9.974118189493335240541442572174e-1
    ],
    // t = -13.75
    [
        F64x2::from_bits(0xBFB0212B746FCD36, 0x3C33571E8A0E76AE), // -6.300613014882547219880242023776e-2
        F64x2::from_bits(0x3FF0F553A6E046DE, 0x3C86FAF483FDDFBE), // 1.059894229755108647668913737521e0
        F64x2::from_bits(0xBFD24FDA7BC15745, 0x3C799A252D3BF3B0), // -2.861238678500408805209077779649e-1
        F64x2::from_bits(0xBFCE928D1D787582, 0xBC5DA3BEDF87DCDC), // -2.388473886393676701589311197890e-1
    ],
    // t = -13.5
    [
        F64x2::from_bits(0x3FC87212C90DA7DE, 0xBC259A72AC7B2272), // 1.909812432962202926851340386850e-1
        F64x2::from_bits(0x3FEA722318F6DF59, 0x3C8DD414ACDB55AC), // 8.264327514252542379390867302024e-1
        F64x2::from_bits(0xBFCCAA65A167EFF5, 0xBC0CD55B63C98221), // -2.239501035800227914673690860743e-1
        F64x2::from_bits(0x3FE652CF94668253, 0x3C688722D914C09F), // 6.976087473340818574176773330646e-1
    ],
    // t = -13.25
    [
        F64x2::from_bits(0x3FD2E6D54238234C, 0x3C7AA0F15E75931C), // 2.953389307663669924904443680813e-1
        F64x2::from_bits(0xBFA88255D561D28C, 0x3C3C133690CD6273), // -4.786937932741350113413401405041e-2
        F64x2::from_bits(0x3F8E108B992FC38C, 0x3C0A328F5A6021BE), // 1.467999516714415410001054455774e-2
        F64x2::from_bits(0x3FF134D57C83D602, 0x3C2700877794D758), // 1.075398908983857371792296378979e0
    ],
    // t = -13
    [
        F64x2::from_bits(0x3FC5F40DD8027025, 0x3C6C89A02ABB84B2), // 1.715104393705370446315817461081e-1
        F64x2::from_bits(0xBFEBE37D3C4A1349, 0xBC88075B1AEB6EFF), // -8.715196778799533667224616334909e-1
        F64x2::from_bits(0x3FCF0DF34643D025, 0x3C4308B1D2FA87CE), // 2.426132290926271993334156622778e-1
        F64x2::from_bits(0x3FE3F069A43D015B, 0xBC61AEFBB6B8C88A), // 6.230972488192877335370886575705e-1
    ],
    // t = -12.75
    [
        F64x2::from_bits(0xBFB54014825DA218, 0x3C14F6C364985EE4), // -8.300903494603806247919307163569e-2
        F64x2::from_bits(0xBFF069709C999192, 0xBC90AB4A4E2C372E), // -1.025742160528441862078425361209e0
        F64x2::from_bits(0x3FD25ABC4CC42927, 0x3C7755974E9FC77C), // 2.867880582495154409338779361654e-1
        F64x2::from_bits(0xBFD29C848E395DE5, 0xBC645F9BCD824152), // -2.908030880839461240299072574829e-1
    ],
    // t = -12.5
    [
        F64x2::from_bits(0xBFD1AE7B7F765332, 0xBC53131BE6869E63), // -2.762745613811602482251711382965e-1
        F64x2::from_bits(0xBFDAD6531395C1A7, 0xBC65C79FBE836B09), // -4.193313304195051644060210937263e-1
        F64x2::from_bits(0x3FBDF5E614E94848, 0x3C508532EE8D2531), // 1.170333672573927766021015095319e-1
        F64x2::from_bits(0xBFEF2F3D4D87F98F, 0xBC81DE1506E70416), // -9.745165361671740721560628257811e-1
    ],
    // t = -12.25
    [
        F64x2::from_bits(0xBFD121173B1A431B, 0x3C715EDFF1180603), // -2.676446988271422982355209581730e-1
        F64x2::from_bits(0x3FDEC698B43E8317, 0x3C384C9ED96E9794), // 4.808713684270044543672564269978e-1
        F64x2::from_bits(0xBFC1C8C7EEE7D95B, 0x3C23908C408718C6), // -1.389398495227379398934660130520e-1
        F64x2::from_bits(0xBFEE11C6CACB3ACF, 0x3C8229EDEE8B48FA), // -9.396699868028351683369815671181e-1
    ],
    // t = -12
    [
        F64x2::from_bits(0xBFB109C28C3CF34F, 0xBC53B7C60D2D38CA), // -6.655517505437312947418966235960e-2
        F64x2::from_bits(0x3FF05EA911169424, 0x3C8A0A3D60604884), // 1.023110453367970729895984322363e0
        F64x2::from_bits(0xBFD2ED1335C9AF37, 0x3C77A043A245C0C5), // -2.957199120780730567294575143568e-1
        F64x2::from_bits(0xBFCE4D3D9BCC24EC, 0xBC41F3959830DBED), // -2.367321978311233163262142861470e-1
    ],
    // t = -11.75
    [
        F64x2::from_bits(0x3FC74C9A0F1CB32F, 0xBC6EF0F7D5C04FAF), // 1.820252012052149958182954436644e-1
        F64x2::from_bits(0x3FEAEE904B2BBEB4, 0xBC8AFEAD1CA7F8D0), // 8.416215389424536817093211127588e-1
        F64x2::from_bits(0xBFCF47A20BA1DAA8, 0xBC61D2C0C3479A31), // -2.443735653684673347487793806398e-1
        F64x2::from_bits(0x3FE3CD54059E9518, 0x3C686C018636B070), // 6.188144788508536215274658057118e-1
    ],
    // t = -11.5
    [
        F64x2::from_bits(0x3FD38C0CC8F233E8, 0x3C6980C2C9C19EB1), // 3.054229700435926563996097722489e-1
        F64x2::from_bits(0x3FB6751715E19E45, 0xBC416D929C1DB8D4), // 8.772415432178444305360574444025e-2
        F64x2::from_bits(0xBF987BAC1B5C1A41, 0x3C0E484132D44EB6), // -2.390927235594575754903595198556e-2
        F64x2::from_bits(0x3FF090B26B8D790E, 0xBC730FD31ED76EE7), // 1.035326404693083440897668357641e0
    ],
    // t = -11.25
    [
        F64x2::from_bits(0x3FCC70B34680CFB0, 0x3C2F1BA11C89B1E4), // 2.221893400434260543268814878728e-1
        F64x2::from_bits(0xBFE6BEB71B7944FA, 0xBC8AEB15CDA49C79), // -7.107806717487215509397269240667e-1
        F64x2::from_bits(0x3FCB4F6F0078022F, 0xBC6AF70DDC179B8A), // 2.133616211213093688167933278352e-1
        F64x2::from_bits(0x3FE8008978808112, 0x3C64598E9A24818E), // 7.500655511430827181986216140299e-1
    ],
    // t = -11
    [
        F64x2::from_bits(0xBF81F08C2B097DEC, 0xBC10E5EFE265218F), // -8.759589255702381289966088468981e-3
        F64x2::from_bits(0xBFF06FEF5A3A8994, 0x3C80319EC6DDC6CF), // -1.027327873664579421461187314031e0
        F64x2::from_bits(0x3FD3D1623AC98142, 0x3C37CE5E03BC1252), // 3.096547674267818863329629932166e-1
        F64x2::from_bits(0xBF968D32328AD716, 0xBC3D94D7A2EF784A), // -2.202299531446446655902909489751e-2
    ],
    // t = -10.75
    [
        F64x2::from_bits(0xBFCE08EDE74B6F30, 0xBC3B0D2B4824D72F), // -2.346475009315951446807426435788e-1
        F64x2::from_bits(0xBFE5AEA832254234, 0x3C7D37A917EC85EE), // -6.775704364209275989126974529540e-1
        F64x2::from_bits(0x3FCA3C663A4C78BE, 0x3C58D9E1A26ADD34), // 2.049682411824899684368642742368e-1
        F64x2::from_bits(0xBFE8783B2BDDD2F0, 0xBC8A81ED52257984), // -7.646766526000572843568955221683e-1
    ],
    // t = -10.5
    [
        F64x2::from_bits(0xBFD3F6989DD42C6F, 0xBC44FD9068658862), // -3.119260350510506008546185721217e-1
        F64x2::from_bits(0x3FB748FD69AB57E9, 0x3C40630F0904C388), // 9.095748739068167287889810875466e-2
        F64x2::from_bits(0xBF9F15ACF29BF878, 0xBC334ECD6608E4A0), // -3.035612326402101317808749158894e-2
        F64x2::from_bits(0xBFF02F923F4C68A3, 0x3C923302DF104792), // -1.011614081630377518647217460742e0
    ],
    // t = -10.25
    [
        F64x2::from_bits(0xBFC902E6C3085CA2, 0xBC52106178EE45D7), // -1.954010441120078195621740064740e-1
        F64x2::from_bits(0x3FE933A16E10DDAA, 0xBC836E6EE2BD8980), // 7.875525617336524782246760079524e-1
        F64x2::from_bits(0xBFCFAC2AD02EB0C3, 0xBC65B2753F289F13), // -2.474416271701383823471756986211e-1
        F64x2::from_bits(0xBFE436F5CCBD0E0B, 0xBC81BD8BAEA5416E), // -6.317090033339253537729648796472e-1
    ],
    // t = -10
    [
        F64x2::from_bits(0x3FA49A7FE67FE71A, 0x3C2CFFAB679ACA86), // 4.024123848644319068943031402993e-2
        F64x2::from_bits(0x3FEFE1673AD658B8, 0x3C6897ABFD0A462D), // 9.962650441327900559045725412889e-1
        F64x2::from_bits(0xBFD423B6DE41FE2B, 0xBC744E432DDC4DDF), // -3.146798296438386331617542115023e-1
        F64x2::from_bits(0x3FBE91EC5FBBC316, 0x3C58447AD9427B63), // 1.194141133999092382775253366815e-1
    ],
    // t = -9.75
    [
        F64x2::from_bits(0x3FD02B010D60B530, 0xBC68CD54F25EF952), // 2.526247625963433555270493147829e-1
        F64x2::from_bits(0x3FE3B70E82101278, 0xBC7669CB69C1FF4C), // 6.160957851685244596069874271940e-1
        F64x2::from_bits(0xBFC8FC6C9FF1E06E, 0xBC6A57C2BE3FE67F), // -1.952033787708872928258887752787e-1
        F64x2::from_bits(0x3FE91624536F990D, 0x3C84B28C3A7E1CDB), // 7.839528684242239945809339484261e-1
    ],
    // t = -9.5
    [
        F64x2::from_bits(0x3FD46C3007401113, 0x3C7868B90FECD717), // 3.191032477191282013757477619471e-1
        F64x2::from_bits(0xBFBBAC228325CE64, 0xBC5D08F5B5F269C8), // -1.080953188118712389963452687619e-1
        F64x2::from_bits(0x3FA3589CB99D6D4C, 0xBC236E38BFC14C1B), // 3.778543248946650226563065894493e-2
        F64x2::from_bits(0x3FEF82C714D23EBE, 0xBC8814C68CEC784B), // 9.847140700021197039206687085971e-1
    ],
    // t = -9.25
    [
        F64x2::from_bits(0x3FCA454C4D3B4444, 0x3C66AA5800C8EE2C), // 2.052398087603554231491397899079e-1
        F64x2::from_bits(0xBFE8295E21B203EF, 0xBC8D7D18176C71C3), // -7.550497682678933243094630425414e-1
        F64x2::from_bits(0x3FD00083AC33E945, 0x3C5688A02292A1C1), // 2.500313932101970220157749538884e-1
        F64x2::from_bits(0x3FE431D8ECA3A41C, 0xBC8D511717380D65), // 6.310848829135724027323436460362e-1
    ],
    // t = -9
    [
        F64x2::from_bits(0xBF96AA38E8BD0844, 0xBC141A9DA60938DE), // -2.213372154734140367416924227415e-2
        F64x2::from_bits(0xBFEF38A3AB3ED723, 0xBC87E5E6BCDEC851), // -9.756639809263315947126596842726e-1
        F64x2::from_bits(0x3FD4CBEFDBCA6EC4, 0x3C681BB74B2B72D3), // 3.249473234552449179194280555745e-1
        F64x2::from_bits(0xBFAD6399A376DCFB, 0xBC3E65A4898337A0), // -5.740051384366925439265490030687e-2
    ],
    // t = -8.75
    [
        F64x2::from_bits(0xBFCE7E52680DF24B, 0xBC6E03BC1B0ADBE7), // -2.382300384596355144189426074851e-1
        F64x2::from_bits(0xBFE5903AD9330ED5, 0xBC8CDFD8103872AE), // -6.738561861206686044625639290644e-1
        F64x2::from_bits(0x3FCCDBB3E9AEEAE0, 0xBC6A1B70D3AD0B1F), // 2.254547968894575646924036338936e-1
        F64x2::from_bits(0xBFE6597F0D12BA8F, 0xBC85662761D334F8), // -6.984248404822483251974646743788e-1
    ],
    // t = -8.5
    [
        F64x2::from_bits(0xBFD52379AA33D405, 0xBC61168FD82B3E62), // -3.302902376302088790217001028989e-1
        F64x2::from_bits(0xBFA08B600C36AC3C, 0x3C3AE639F9B56D17), // -3.231334828463913587288273852933e-2
        F64x2::from_bits(0x3F7FC31DAF10CE95, 0x3C0CC53868499B00), // 7.754436447658404431949082773798e-3
        F64x2::from_bits(0xBFEED0A4B5E7EC22, 0x3C8293F17BAB47C5), // -9.629691651201747981359277936848e-1
    ],
    // t = -8.25
    [
        F64x2::from_bits(0xBFD04A52B594C2BE, 0x3C4EF31C32E06CA4), // -2.545363209965606465541003012341e-1
        F64x2::from_bits(0x3FE378FB5D04C51C, 0x3C785543EE467EC5), // 6.085182968874138997985681197191e-1
        F64x2::from_bits(0xBFCB74190C928B3E, 0x3C178897D0BC3157), // -2.144805251492360450712715481479e-1
        F64x2::from_bits(0xBFE79BFB812CF86C, 0x3C89B05E4E34E99A), // -7.377908251726358999767283817116e-1
    ],
    // t = -8
    [
        F64x2::from_bits(0xBFAAFC28073ABB84, 0xBC3291D25865C5A0), // -5.270505035638620262208267579389e-2
        F64x2::from_bits(0x3FEDF01D7E1F41FA, 0x3C551542FAB32A15), // 9.355609381983065510255224621326e-1
        F64x2::from_bits(0xBFD53339D484C3A1, 0x3C52CAE57E4BE869), // -3.312515807511378599698762392757e-1
        F64x2::from_bits(0xBFC468DFB8B805B9, 0x3C6AE4F22C93CDE6), // -1.594504978129813893499357336502e-1
    ],
    // t = -7.75
    [
        F64x2::from_bits(0x3FC665AD04B8F079, 0x3C3ECDE40CD7CF59), // 1.749779007967651473000206903001e-1
        F64x2::from_bits(0x3FE9F59E5A8E3EBA, 0xBC8AE51B7785C085), // 8.112327355065282552278205899944e-1
        F64x2::from_bits(0xBFD2839ED9101846, 0x3C6646F9BECDF32E), // -2.892834777597993359798399275399e-1
        F64x2::from_bits(0x3FDE9702D3352B1A, 0xBC72E3E40536B2E7), // 4.779669821333968174248782578633e-1
    ],
    // t = -7.5
    [
        F64x2::from_bits(0x3FD497F92CA01E61, 0xBC72B7704D5BE8A0), // 3.217757163806478752673285436798e-1
        F64x2::from_bits(0x3FD4675FFD3B25D9, 0xBC7F4055FB45E268), // 3.188095066985545962100629060794e-1
        F64x2::from_bits(0xBFBCCA682E65D075, 0xBC5578B3D6BCF263), // -1.124634850764908063843208150544e-1
        F64x2::from_bits(0x3FEC16F4CF72FBB8, 0x3C7B85BC4F9B20DB), // 8.778022815457609223675812589159e-1
    ],
    // t = -7.25
    [
        F64x2::from_bits(0x3FD4B82A61952443, 0xBC65A02FBD23D475), // 3.237405732111861462212963107673e-1
        F64x2::from_bits(0xBFD336F3ACC28469, 0xBC74929F73D76E71), // -3.002289950473540814628952839942e-1
        F64x2::from_bits(0x3FBD97638DBD4C4B, 0x3C3A8AD61B8A67C8), // 1.155912610095565660151131034532e-1
        F64x2::from_bits(0x3FEC086D5EAE5CEA, 0xBC8928AB272EA384), // 8.760287141075455260455826158286e-1
    ],
    // t = -7
    [
        F64x2::from_bits(0x3FC79683B0571A28, 0xBC6899B4D6652337), // 1.842808352505056372799415198167e-1
        F64x2::from_bits(0xBFE8AC195288A6BD, 0xBC676C5A4DF35A4A), // -7.710081684101265477312516545347e-1
        F64x2::from_bits(0x3FD2CCFF6EDADFCF, 0x3C7A16865856A333), // 2.937620718544140201236471886909e-1
        F64x2::from_bits(0x3FDFE33D46ED5876, 0xBC6CBBE9411D379A), // 4.982445900581134887461169261108e-1
    ],
    // t = -6.75
    [
        F64x2::from_bits(0xBFA117CFAF866D7D, 0xBC33C44458491249), // -3.338479058876495899085204876535e-2
        F64x2::from_bits(0xBFED03B837256366, 0xBC84858361F7222B), // -9.067040516921281235350748964295e-1
        F64x2::from_bits(0x3FD64B3805C6158A, 0xBC464ED64501B5B0), // 3.483409935364184504508073545535e-1
        F64x2::from_bits(0xBFB2EC35A8E48A1C, 0x3C51002B7E690B20), // -7.391677258832668093675208709642e-2
    ],
    // t = -6.5
    [
        F64x2::from_bits(0xBFCE7773026E4ABD, 0xBC60DF255161DFB7), // -2.380203019971158035944441034961e-1
        F64x2::from_bits(0xBFE59935F836551A, 0x3C89E98FA3254A79), // -6.749524925132021729989387543664e-1
        F64x2::from_bits(0x3FD0B46E6F1901F0, 0xBC5C0EB7B7ECBB88), // 2.610126576364839518174202614884e-1
        F64x2::from_bits(0xBFE31C05A83B4C27, 0xBC8DB41137FE325B), // -5.971706662916220169762741441342e-1
    ],
    // t = -6.25
    [
        F64x2::from_bits(0xBFD6600B39FD3304, 0xBC7DB1B4E70A0C12), // -3.496120516108905098546429475489e-1
        F64x2::from_bits(0xBFC87583B9AF1C1B, 0xBC4CA5C1126710E8), // -1.910862595234171543685577403638e-1
        F64x2::from_bits(0x3FB2210E6A9DA593, 0x3C406F3352F3F6CC), // 7.081689932751649053090894183665e-2
        F64x2::from_bits(0xBFEBE574E9E20D20, 0x3C87FC89F6FDB458), // -8.717598503139107474448063207871e-1
    ],
    // t = -6
    [
        F64x2::from_bits(0xBFD510B6EB1815D2, 0xBC660A0EEAEB9601), // -3.291451736298231052314485825290e-1
        F64x2::from_bits(0x3FD623CE99198C5B, 0x3C71E8F0D2DB460F), // 3.459354872813428949297794348338e-1
        F64x2::from_bits(0xBFC2C7032D16920B, 0x3C6237848598ED80), // -1.466983766705570378752607399885e-1
        F64x2::from_bits(0xBFEA0344501F8029, 0xBC82D7BDAA2B78A0), // -8.128987851050670004246809994469e-1
    ],
    // t = -5.75
    [
        F64x2::from_bits(0xBFC82BFA57A7C26C, 0x3C57E958EBF16EAD), // -1.888420989994473668025316562222e-1
        F64x2::from_bits(0x3FE7A73ECC8B9B11, 0x3C8C245F7D4BD541), // 7.391656870866844463963154600278e-1
        F64x2::from_bits(0xBFD3EE22630FD44E, 0x3C596ADA807FE807), // -3.114095656777110477216482868990e-1
        F64x2::from_bits(0xBFDDDDE4B3A404DA, 0xBC4503FAEBD42190), // -4.666682962707234839120346691220e-1
    ],
    // t = -5.5
    [
        F64x2::from_bits(0x3F92355309057E0A, 0xBC236A0987288B82), // 1.778154127657497560302015149724e-2
        F64x2::from_bits(0x3FEBA780EC73EA42, 0x3C86927A0BB3EC1B), // 8.641972177713983907721118946923e-1
        F64x2::from_bits(0xBFD78A4170E631BF, 0x3C7FC1613889F9E9), // -3.678134539157119910947077770812e-1
        F64x2::from_bits(0x3F99B6D9D03547E0, 0xBC3D13C3A7DCED1B), // 2.511158307363092598875579956098e-2
    ],
    // t = -5.25
    [
        F64x2::from_bits(0x3FCC08806801325D, 0xBC6CC5F22762F693), // 2.190094478450132095664432197169e-1
        F64x2::from_bits(0x3FE6733C101C4B90, 0xBC8DC74691D81F28), // 7.015667261751889521539581082369e-1
        F64x2::from_bits(0xBFD34945F2F5526C, 0x3C645B79E6BC086A), // -3.013472435607471489581659430925e-1
        F64x2::from_bits(0x3FDF3CBE8A417766, 0x3C5912CF144F63B3), // 4.880825376657099755574737483760e-1
    ],
    // t = -5
    [
        F64x2::from_bits(0x3FD672DE4D9E1D32, 0xBC70897D7849497D), // 3.507610090241143197880163276967e-1
        F64x2::from_bits(0x3FD4F0BA25CB5A72, 0xBC7455D5948EC1AC), // 3.271928185544431367948786774266e-1
        F64x2::from_bits(0xBFC1B6146E96CED2, 0xBC03F7020F00518D), // -1.383691349016005768500291756026e-1
        F64x2::from_bits(0x3FE8E8BFCE7BAA41, 0x3C7031517AD37FB2), // 7.784117730018992460944232099040e-1
    ],
    // t = -4.75
    [
        F64x2::from_bits(0x3FD80F453C3D9834, 0x3C539D5533E21C10), // 3.759320343291421327236144976706e-1
        F64x2::from_bits(0xBFC044CCC5FEDAD9, 0x3C5E9DAD6C93324D), // -1.270996062064202669853759238881e-1
        F64x2::from_bits(0x3FB135B410B250B1, 0x3C2575BF7614F184), // 6.722569854383909973441449601051e-2
        F64x2::from_bits(0x3FEA5E27782FB5F1, 0xBC80F541BD63AB55), // 8.239934298887288938207610487549e-1
    ],
    // t = -4.5
    [
        F64x2::from_bits(0x3FD2B2A1940487E5, 0xBC7790E99732600B), // 2.921527810559594668815688954853e-1
        F64x2::from_bits(0xBFE0BF62C807EEA1, 0xBC504AA4C9D53661), // -5.233625323157477007084954792739e-1
        F64x2::from_bits(0x3FD03F731A8A0B1D, 0x3C6CD6A31964F4B1), // 2.538726576969326368005244553070e-1
        F64x2::from_bits(0x3FE44FD4425CA3D8, 0xBC795429F09D3B8A), // 6.347447677736637097333254122672e-1
    ],
    // t = -4.25
    [
        F64x2::from_bits(0x3FC05B30E2B75DFE, 0xBC0B3784E9606088), // 1.277829272282672843739985663937e-1
        F64x2::from_bits(0xBFE84BEB2BF68F2D, 0xBC86F303623569E5), // -7.592674120573740646580717870614e-1
        F64x2::from_bits(0x3FD7C16237E4164B, 0x3C638D38B28516DA), // 3.711782022295195356124806943330e-1
        F64x2::from_bits(0x3FD246307DAF268B, 0x3C784C8443A699A5), // 2.855340220818127427809932712489e-1
    ],
    // t = -4
    [
        F64x2::from_bits(0xBFB1FCEC060D9F19, 0xBC3AAB18CC944A6D), // -7.026553294928951509908431163180e-2
        F64x2::from_bits(0xBFE94CD44C4FB752, 0xBC8BB5AEC3FBA95A), // -7.906285753685813802964544458279e-1
        F64x2::from_bits(0x3FD91A5F98578929, 0xBC76B828106F918E), // 3.922347057069992895544918276468e-1
        F64x2::from_bits(0xBFBDDE1F4F8DAE06, 0x3C21CB7910FC4AA2), // -1.166705674383408936795672429767e-1
    ],
    // t = -3.75
    [
        F64x2::from_bits(0xBFD01A6C2AAD1263, 0xBC4CBF8B51C7F3E6), // -2.516127030142227303326977416092e-1
        F64x2::from_bits(0xBFE43D1019AA26D4, 0xBC72FA63403F7DF2), // -6.324539662611763533324433970960e-1
        F64x2::from_bits(0x3FD44CC41D6706C7, 0xBC450522F65F2F56), // 3.171854292996666967977096344660e-1
        F64x2::from_bits(0xBFDDF074180BBD9C, 0x3C7C41A913403E99), // -4.678011164496298256716635508031e-1
    ],
    // t = -3.5
    [
        F64x2::from_bits(0xBFD808BF043B852C, 0x3C77ECACA90DFC33), // -3.755338231404319119343969515802e-1
        F64x2::from_bits(0xBFD5FAFA2AAD827C, 0x3C5D193E9A4E3BC1), // -3.434434334540481462879373740987e-1
        F64x2::from_bits(0x3FC59FD2127C817E, 0xBC6070A30979F499), // 1.689398374810586118434427695409e-1
        F64x2::from_bits(0xBFE62E0234000EF1, 0x3C7C9BDF0CD48A3F), // -6.931162849072888017524436126706e-1
    ],
    // t = -3.25
    [
        F64x2::from_bits(0xBFDAD11D0561308F, 0x3C768D00F6D9AFD0), // -4.190132668052308022390469206240e-1
        F64x2::from_bits(0xBF641A17B6EE4670, 0xBC09FD6192B2C183), // -2.453848187948186497374622876146e-3
        F64x2::from_bits(0xBF906B1AFC915402, 0xBC1040FC16AB75E2), // -1.603357473898726231854718112212e-2
        F64x2::from_bits(0xBFE84FF2C13DDB42, 0xBC821692F48FF851), // -7.597593092203640855051072562928e-1
    ],
    // t = -3
    [
        F64x2::from_bits(0xBFD83E7E4EA6959A, 0xBC65B74FD6B2805F), // -3.788142936776580743472439164997e-1
        F64x2::from_bits(0x3FD42223F628D022, 0x3C5C2FC6A8EB3ADA), // 3.145837692165988136507872660659e-1
        F64x2::from_bits(0xBFC9618DF2354B90, 0xBC6C797A42379A6C), // -1.982896263749265432206448545725e-1
        F64x2::from_bits(0xBFE59E9B6D47CCDF, 0xBC66A16C5E0508E7), // -6.756112226852585376680320451823e-1
    ],
    // t = -2.75
    [
        F64x2::from_bits(0xBFD12EF2F87EB4E5, 0xBC47091340E0D3CF), // -2.684905459125970808571045227290e-1
        F64x2::from_bits(0x3FE1A48FBEBFF8B1, 0x3C874FA5C2D217FE), // 5.513380742629775803904987719004e-1
        F64x2::from_bits(0xBFD60A4111CA8EBD, 0xBC58C5ABDD979DF2), // -3.443758653395255274069067482088e-1
        F64x2::from_bits(0xBFDE9DE415E3B6DA, 0x3C6B138C7CDDEA7F), // -4.783868993534788776315024621769e-1
    ],
    // t = -2.5
    [
        F64x2::from_bits(0xBFBCC155EC43247D, 0xBC36CEC3F90EA51B), // -1.123250676929660891874631001402e-1
        F64x2::from_bits(0x3FE5B9295E8EF584, 0x3C8E4DC78C6E4C30), // 6.788527342647943633721400308225e-1
        F64x2::from_bits(0xBFDBACCF4DA71FB5, 0x3C77BFFBDA2CB773), // -4.324224718407052930284195036923e-1
        F64x2::from_bits(0xBFCC36BA46486BC7, 0x3C67C7B2C537089A), // -2.204201548746295876833984275337e-1
    ],
    // t = -2.25
    [
        F64x2::from_bits(0x3FAF89DC01D9726A, 0xBC4832218DB0E43F), // 6.159865877700527751717642376258e-2
        F64x2::from_bits(0x3FE63D92A0BF2328, 0x3C7469A2A92B421D), // 6.950162067015286559394633840143e-1
        F64x2::from_bits(0xBFDD0D095A247490, 0x3C6C2C4714262C3E), // -4.539206867501173070455352154418e-1
        F64x2::from_bits(0x3FA780C9A0C98F7A, 0x3C3D1707DE4439AC), // 4.590444648491050374562271090557e-2
    ],
    // t = -2
    [
        F64x2::from_bits(0x3FCD1BAFC57F31D0, 0xBC5112970239D372), // 2.274074282016855759919244360379e-1
        F64x2::from_bits(0x3FE3C8C724515C8F, 0x3C602AA022B440DC), // 6.182590207416910414062642913325e-1
        F64x2::from_bits(0xBFDA632A64D50C9E, 0x3C7DFF580C98A00B), // -4.123025879563984880832340546115e-1
        F64x2::from_bits(0x3FD1D7C7AF0D64E9, 0x3C670153B5258AAA), // 2.787951669211695226850975694110e-1
    ],
    // t = -1.75
    [
        F64x2::from_bits(0x3FD76413DDDFBBE9, 0x3C659EE86BFE6F60), // 3.654832522142315669961822981181e-1
        F64x2::from_bits(0x3FDEA23A339F5D7C, 0xBC7E8E89D31C27C3), // 4.786515716673062927094621707797e-1
        F64x2::from_bits(0xBFD472BD4AD1A9A3, 0xBC7A938F4C9CF1C6), // -3.195031386045689933559551075291e-1
        F64x2::from_bits(0x3FDCF5AC02C6404B, 0x3C4F0F252F929A1E), // 4.524946238607340831906724688596e-1
    ],
    // t = -1.5
    [
        F64x2::from_bits(0x3FDDB661389897F8, 0xBC62960EFDD05A66), // 4.642565777488694064742733669192e-1
        F64x2::from_bits(0x3FD3C9B8221EF635, 0x3C766599A5A064E2), // 3.091869672024104204161689166457e-1
        F64x2::from_bits(0xBFC88C680544C949, 0xBC5A6984E3C43758), // -1.917848611570412200129417232656e-1
        F64x2::from_bits(0x3FE1DA621814E0DA, 0x3C8BD4F9A128F4D9), // 5.579081030218973541316939627343e-1
    ],
    // t = -1.25
    [
        F64x2::from_bits(0x3FE0A43669C08DDA, 0x3C75C474F0CEB0E9), // 5.200454774352991826955478732193e-1
        F64x2::from_bits(0x3FC1CD5BF011B8A1, 0xBC6A3DAD35DFBA26), // 1.390795633519177507127373116956e-1
        F64x2::from_bits(0xBFA77BF0DC44EA76, 0xBC37135D81C5C682), // -4.586746872742690620177231141081e-2
        F64x2::from_bits(0x3FE331AD892034B6, 0xBC3B64154316ED42), // 5.998141935575833623845457522842e-1
    ],
    // t = -1
    [
        F64x2::from_bits(0x3FE1235093D83DA5, 0x3C897F341F179F09), // 5.355608832923521187995165656389e-1
        F64x2::from_bits(0xBF84CF103BCC6624, 0x3C14862738CF8E58), // -1.016056711664520939504546984536e-2
        F64x2::from_bits(0x3FBA9F92AAC23D81, 0x3C5AAFA2DF2B95AC), // 1.039973894969446118886899909786e-1
        F64x2::from_bits(0x3FE2F4BDBACDDE78, 0xBC8CBBA0D7B5A7A6), // 5.923756264227923508167792291816e-1
    ],
    // t = -0.75
    [
        F64x2::from_bits(0x3FE09197CAA1560C, 0x3C7703BE9677A74E), // 5.177725751515836116525029549064e-1
        F64x2::from_bits(0xBFC020755035F8BB, 0x3C67333A5C7BC687), // -1.259905473379541906292026393536e-1
        F64x2::from_bits(0x3FCFB73F063E539D, 0x3C6FC13F04FE742B), // 2.477797298894558879525598698501e-1
        F64x2::from_bits(0x3FE1BE427DF752A7, 0xBC8D4CA39EE4E306), // 5.544750652575956608944316344972e-1
    ],
    // t = -0.5
    [
        F64x2::from_bits(0x3FDE72543CD05E5D, 0x3C5B0BA76C6604EB), // 4.757280916105395887986437782813e-1
        F64x2::from_bits(0xBFCA1F5921E923F0, 0xBC4CED6F2340462C), // -2.040816703395473861448172017949e-1
        F64x2::from_bits(0x3FD857B2AEA4F3FD, 0xBC7418073D7535CB), // 3.803526597510538501697124937263e-1
        F64x2::from_bits(0x3FE0309BE63EAA03, 0xBC85242C1404908D), // 5.059337136238471665702604378969e-1
    ],
    // t = -0.25
    [
        F64x2::from_bits(0x3FDACC625315F862, 0xBC64BC34FB4EF37B), // 4.187246142754529242283811576924e-1
        F64x2::from_bits(0xBFCF89AE5465C577, 0xBC3974E9E9CF8049), // -2.463891899201759730286849771583e-1
        F64x2::from_bits(0x3FE00B77BF589ACB, 0xBC853DD305987733), // 5.013998734692333889674857802277e-1
        F64x2::from_bits(0x3FDDC50ABF8600D4, 0x3C61777270B5C765), // 4.651514883371537032706729631105e-1
    ],
    // t = 0
    [
        F64x2::from_bits(0x3FD6B8C7962715B8, 0x3C77A96D7BB04E65), // 3.550280538878172392600631860042e-1
        F64x2::from_bits(0xBFD0907F42B70F8B, 0x3C7D1459035AFDE2), // -2.588194037928067984051835601892e-1
        F64x2::from_bits(0x3FE3AD7A9B4A3EA9, 0x3C8D5765B40267BD), // 6.149266274460007351509223690936e-1
        F64x2::from_bits(0x3FDCB0C1A680C8A1, 0xBC7D3DE8103B7766), // 4.482883573538263579148237103988e-1
    ],
    // t = 0.25
    [
        F64x2::from_bits(0x3FD2A26E236CDDD1, 0x3C6949A6057C7D76), // 2.911639543485452062721071941866e-1
        F64x2::from_bits(0xBFCFE1446CDDEA9C, 0xBC399993192B3AD0), // -2.490621120048971418037222001683e-1
        F64x2::from_bits(0x3FE751E506EEF419, 0x3C6BB9FF884C2997), // 7.287469039362150078694543181809e-1
        F64x2::from_bits(0x3FDE1234AF392559, 0xBC73102E8B7B90B9), // 4.698611937679593565451907822836e-1
    ],
    // t = 0.5
    [
        F64x2::from_bits(0x3FCDA822D7438440, 0x3C654C9A822E0B1E), // 2.316936064808334897691252545099e-1
        F64x2::from_bits(0xBFCCC9DE4B290E91, 0xBC5982172F63597C), // -2.249105326646838931359969903286e-1
        F64x2::from_bits(0x3FEB563CCF3B4098, 0x3C8D609ADA7A2C22), // 8.542770431031554933000487987952e-1
        F64x2::from_bits(0x3FE16D2371290BEE, 0x3C5756B6731BBFAD), // 5.445725641405923018271640182178e-1
    ],
    // t = 0.75
    [
        F64x2::from_bits(0x3FC6F47DF7821461, 0xBC45B112FB3A1451), // 1.793363054786452336149763775180e-1
        F64x2::from_bits(0xBFC8B9F7189A67BF, 0x3C6854C627E188F1), // -1.931752081043764562813758971023e-1
        F64x2::from_bits(0x3FF01C63959E9464, 0x3C87229CCA4F6B30), // 1.006930908633216366242233100406e0
        F64x2::from_bits(0x3FE616EF66F630BA, 0xBC703C4E8680CEB8), // 6.902997027368862196405423897753e-1
    ],
    // t = 1
    [
        F64x2::from_bits(0x3FC151430BBAF656, 0x3C3DFFCFECF529DA), // 1.352924163128814155241474235155e-1
        F64x2::from_bits(0xBFC45EF17FCE4FAE, 0xBC69819033DAAA41), // -1.591474412967932127875002524972e-1
        F64x2::from_bits(0x3FF3519B674BDC80, 0xBC9346A88B7545F4), // 1.207423594952871259436378817028e0
        F64x2::from_bits(0x3FEDD683E1F130CA, 0xBC600A986599CA17), // 9.324359333927756329594514536744e-1
    ],
    // t = 1.25
    [
        F64x2::from_bits(0x3FB9824E0CF4AD0A, 0x3C3F096E49004E1A), // 9.964454475691667147460101896242e-2
        F64x2::from_bits(0xBFC030B6AD9CD6C7, 0x3C31FDDA04806F5F), // -1.264866206853893772182923158977e-1
        F64x2::from_bits(0x3FF7C00DEB9D13BD, 0x3C761131C06B5128), // 1.484388275495106121182301537941e0
        F64x2::from_bits(0x3FF4F697ECF37959, 0xBC8ABDFF3408957F), // 1.310203481283300950516580471014e0
    ],
    // t = 1.5
    [
        F64x2::from_bits(0x3FB25E2CCF277DC1, 0xBC566621DB196BA7), // 7.174949700810540967355541648968e-2
        F64x2::from_bits(0xBFB8EE0710605791, 0xBC4E47122ED332CB), // -9.738201284230131921848421820245e-2
        F64x2::from_bits(0x3FFE1024F75B1BBB, 0x3C95BA6AE67E64E8), // 1.878941503747895000909335049476e0
        F64x2::from_bits(0x3FFE2DECE6BE3255, 0x3C80E80AB11ED8B3), // 1.886212254848165488692347024025e0
    ],
    // t = 1.75
    [
        F64x2::from_bits(0x3FA9E44BA0726F86, 0xBC36C8F433F22598), // 5.056988080579487164484523641890e-2
        F64x2::from_bits(0xBFB2A68A7E0003BD, 0x3BF8105B788192EA), // -7.285371376202838525522367052363e-2
        F64x2::from_bits(0x40039E40190E65F5, 0xBCA68606AD0CD075), // 2.452270694496059565667821449609e0
        F64x2::from_bits(0x400617B82989272C, 0x3C6206403C0CC11E), // 2.761581730363920080862692053922e0
    ],
    // t = 2
    [
        F64x2::from_bits(0x3FA1E1935C04B623, 0x3C277E337C777CA1), // 3.492413042327437913532208079181e-2
        F64x2::from_bits(0xBFAB2EA9B1BFCCDC, 0xBC20E11773C85895), // -5.309038443365363170399918587870e-2
        F64x2::from_bits(0x400A627FA1A00F78, 0xBC8BD86E27529F2F), // 3.298094999978214710280604425223e0
        F64x2::from_bits(0x4010671931FF0625, 0x3CA27BE74A23118F), // 4.100682049932889889382034079178e0
    ],
    // t = 2.25
    [
        F64x2::from_bits(0x3F9838ED439D41EE, 0xBC373F326BB61937), // 2.365465855774744620687506627587e-2
        F64x2::from_bits(0xBFA35517671BDC20, 0x3C339D97D90499B0), // -3.775857099201851312631718748093e-2
        F64x2::from_bits(0x401240B907888320, 0xBCB5870643B99A72), // 4.563205831248324940604419737992e0
        F64x2::from_bits(0x4018B0B313AEA784, 0x3C9F25624D51ED86), // 6.172558124098831308458720859432e0
    ],
    // t = 2.5
    [
        F64x2::from_bits(0x3F901A74DA795DF8, 0xBC3661768DCC2996), // 1.572592338047048999526604654076e-2
        F64x2::from_bits(0xBF9AE182CE295C41, 0x3C3E011DD594AE4A), // -2.625088103590323036489549629723e-2
        F64x2::from_bits(0x4019ED3878FDD253, 0xBC8FF0613D3B03AF), // 6.481660738460578608072612957487e0
        F64x2::from_bits(0x4022D7C4CC0B7200, 0xBC94DB22BA6C82D5), // 9.421423317334301755582308885728e0
    ],
    // t = 2.75
    [
        F64x2::from_bits(0x3F85080603DD42B5, 0xBC1D8022FEA77498), // 1.026920985501198752263947149980e-2
        F64x2::from_bits(0xBF924AF70A194640, 0x3C276C3F4856CD28), // -1.786409377229447529095519052370e-2
        F64x2::from_bits(0x4022DD60C8A449A8, 0x3CC9DA6AE808A8AB), // 9.432379026465086020813815057947e0
        F64x2::from_bits(0x402D2D24AA20366D, 0x3CCE542D932A8117), // 1.458817035334797202830227117706e1
    ],
    // t = 3
    [
        F64x2::from_bits(0x3F7AFF4F7FBD1F1B, 0xBBF36861BAA1ECCC), // 6.591139357460719144257448407961e-3
        F64x2::from_bits(0xBF8865D4AACCF1E1, 0x3C10FA5CC8BBEB11), // -1.191297670595131847376323259302e-2
        F64x2::from_bits(0x402C131CC82CD470, 0xBC90BC2B1834DC62), // 1.403732896373023203174026731378e1
        F64x2::from_bits(0x4036EC1647B07E40, 0xBCCEDEAC4E28A06F), // 2.292221496638217018510473104331e1
    ],
    // t = 3.25
    [
        F64x2::from_bits(0x3F710A8D885744BC, 0x3C10857FF250A8AD), // 4.160454618117256449714454041763e-3
        F64x2::from_bits(0xBF7FEB39BCB22D41, 0xBC0844A9ED151C33), // -7.792687926790721119475965283362e-3
        F64x2::from_bits(0x403554B62FD58B97, 0xBCDE5B0B17BB293B), // 2.133090495074756107986912949873e1
        F64x2::from_bits(0x404247055FAB38A3, 0x3CD81F8A852EC6AF), // 3.655485149250423538420703424286e1
    ],
    // t = 3.5
    [
        F64x2::from_bits(0x3F652B3F78F3BE24, 0x3BFBDD186EA2388A), // 2.584098786989634963277144783300e-3
        F64x2::from_bits(0xBF747F82253F7EF5, 0xBBF51E443B7A98CE), // -5.004413967952582832030249678838e-3
        F64x2::from_bits(0x4040871AD867E1B7, 0x3CD10089598C7541), // 3.305550675461147941425731298192e1
        F64x2::from_bits(0x404D95086C8E04D8, 0x3CE237BE9E5276A1), // 5.916431958136098703457421217952e1
    ],
    // t = 3.75
    [
        F64x2::from_bits(0x3F59E34D16A4BEF0, 0x3BFA31D5B3593218), // 1.580071717921013257845719138246e-3
        F64x2::from_bits(0xBF69DDC9D341F712, 0xBC024D7CA5CB0862), // -3.157514753239784192030096292774e-3
        F64x2::from_bits(0x404A17745BCAA456, 0x3CEBD45A1D2AB3B9), // 5.218323848147035757588652388124e1
        F64x2::from_bits(0x40584B14D5CD7609, 0xBCFFAF388AF17F6A), // 9.717314667763287879055222354540e1
    ],
    // t = 4
    [
        F64x2::from_bits(0x3F4F2E4BCF7C4970, 0x3BE1455B0386DC1C), // 9.515638512048018736214999689001e-4
        F64x2::from_bits(0xBF600B915A6C6845, 0x3BFE1F643915ABF9), // -1.958640950204178900138140918409e-3
        F64x2::from_bits(0x4054F6366AFF2FD4, 0x3CD0A93C4EBE02CB), // 8.384707140846813992258049046105e1
        F64x2::from_bits(0x40643DA7642A41D5, 0xBCFA2604D6D26135), // 1.619266835046134018430949242855e2
    ],
    // t = 4.25
    [
        F64x2::from_bits(0x3F42808AD03B7838, 0xBBDA89A621A0DB8E), // 5.646398353425013377819267939638e-4
        F64x2::from_bits(0xBF53950DBDC09F1E, 0xBBF3B56D41E6F28F), // -1.195205134544914304407708132072e-3
        F64x2::from_bits(0x406120AEDDC9D5CF, 0x3CF1413CD0BFA31B), // 1.370213459913343039830633739948e2
        F64x2::from_bits(0x40711B2E767EA9DB, 0xBD0AA61C214CC656), // 2.736988434741776240915515168506e2
    ],
    // t = 4.5
    [
        F64x2::from_bits(0x3F35A4AE56C7E071, 0x3BDB3F2F3FDE7517), // 3.302503235143089836587325900993e-4
        F64x2::from_bits(0xBF4785E6B71C4DDF, 0xBBE1FDD52F14E1AD), // -7.178665675575088886935542984668e-4
        F64x2::from_bits(0x406C72D190FF68B8, 0xBCF6A5EB5C4A6721), // 2.275880818355997184614108860545e2
        F64x2::from_bits(0x407D522946D820BF, 0x3D1E05C2858C8AA6), // 4.691350773279663979509196771451e2
    ],
    // t = 4.75
    [
        F64x2::from_bits(0x3F28F6D37A3790EA, 0x3BB7876DF6483221), // 1.904614592681605127238217147551e-4
        F64x2::from_bits(0xBF3BD37BB7204CF5, 0x3BD8458B31811727), // -4.245926894565620827979542676971e-4
        F64x2::from_bits(0x4077FFE390F1039E, 0xBD11923ACCAC6D6F), // 3.839930581488241458288171254946e2
        F64x2::from_bits(0x408979D0007370AC, 0x3D1959E50B7AF435), // 8.152265633600959743098279348440e2
    ],
    // t = 5
    [
        F64x2::from_bits(0x3F1C66DF1A2952D5, 0xBBB717141F398254), // 1.083444281360744173498650250335e-4
        F64x2::from_bits(0xBF3036EA91E217E0, 0x3BB27751D4F01515), // -2.474138908684624760002361720631e-4
        F64x2::from_bits(0x40848E561B412157, 0x3D26D2F201663E58), // 6.577920441711711824410805788744e2
        F64x2::from_bits(0x40966F46BCF6F424, 0x3D347BA17D1BFBED), // 1.435819080217982518671721238005e3
    ],
    // t = 5.25
    [
        F64x2::from_bits(0x3F0FE1CBA007C5AD, 0x3BA18B28A6BDF291), // 6.081011452242365287333986492002e-5
        F64x2::from_bits(0xBF229FE777B49F0B, 0xBBCC99D49CED7CE1), // -1.420946171972681576101851025903e-4
        F64x2::from_bits(0x4091DE1B0CD3CC83, 0xBD1A7FF1DE9830F3), // 1.143526416119916324643149265023e3
        F64x2::from_bits(0x40A404D61096B940, 0xBD4D2B207B2D8292), // 2.562418095312226589268331812708e3
    ],
    // t = 5.5
    [
        F64x2::from_bits(0x3F01A92A8107B6EB, 0x3BA2EEE3726B0BF9), // 3.368531190859981442528973405694e-5
        F64x2::from_bits(0xBF1517CE89672D6A, 0x3B9ED6B7CFA207E5), // -8.046339130556514337967075505770e-5
        F64x2::from_bits(0x409F8251F5A78469, 0x3D3F600A548AEE67), // 2.016580038659531394441009506935e3
        F64x2::from_bits(0x40B2188DC1747AE6, 0x3D2779D2EC43FE36), // 4.632553733139042420454022913799e3
    ],
    // t = 5.75
    [
        F64x2::from_bits(0x3EF350EA6647CC50, 0x3B97B565986A6895), // 1.842124619773024582063210167374e-5
        F64x2::from_bits(0xBF078FD3FCC0D042, 0x3B9B9171243A21D7), // -4.494062122298348062874345443598e-5
        F64x2::from_bits(0x40AC2C178113BBC8, 0x3D453814F3EBB68D), // 3.606045906654999422038583493525e3
        F64x2::from_bits(0x40C0911460C99F12, 0x3D4E22FD309D8811), // 8.482159203722640300332186820739e3
    ],
    // t = 6
    [
        F64x2::from_bits(0x3EE4DCA0B3CC0F9D, 0x3B8C4D8338215DB8), // 9.947694360252889570238847668829e-6
        F64x2::from_bits(0xBEF9F7DB9CCFD7A0, 0x3B8BF5AE5C1C3943), // -2.476520039703495475418182538699e-5
        F64x2::from_bits(0x40B9887233ECC0F4, 0xBD552EBDC699901D), // 6.536446104809863453758350024619e3
        F64x2::from_bits(0x40CEB6CD22B725A0, 0x3D3AEF5D38DECF91), // 1.572560262193047683942032295732e4
    ],
    // t = 6.25
    [
        F64x2::from_bits(0x3ED641202C0A3CAA, 0xBB7F58ED27BBD934), // 5.305861748752081026322708937454e-6
        F64x2::from_bits(0xBEEC3F2CDB61EE61, 0xBB6026FCBCB15C36), // -1.346911345145098343914944543104e-5
        F64x2::from_bits(0x40C7731C70F7652E, 0xBD69B5011C84D824), // 1.200622219746055925193535902482e4
        F64x2::from_bits(0x40DCD27A2222CF7E, 0x3D63CFCF89FE88C0), // 2.951390833349478691081232410726e4
    ],
    // t = 6.5
    [
        F64x2::from_bits(0x3EC7741C92B83C35, 0x3B324E59C74E1A32), // 2.795882343204913585459995748811e-6
        F64x2::from_bits(0xBEDE553A2F48A090, 0x3B7C9E8CC569093D), // -7.231931466601792559814248837755e-6
        F64x2::from_bits(0x40D5D126E4DBB412, 0x3D6BB68D4ACF9DA5), // 2.234060771839699815794499069002e4
        F64x2::from_bits(0x40EB5FCFDDF12375, 0xBD6A3F5A32D34BE7), // 5.606249584252286074821910033153e4
    ],
    // t = 6.75
    [
        F64x2::from_bits(0x3EB86CAB0A601B37, 0x3B48E0ACBFF89DDF), // 1.455812744578875868998232086113e-6
        F64x2::from_bits(0xBED015372B187171, 0xBB57BE5609AE4748), // -3.834455740949934238658678740430e-6
        F64x2::from_bits(0x40E48E8C24C1592F, 0xBD806544AA97E090), // 4.210037948672693705753101587476e4
        F64x2::from_bits(0x40FA4EFA1926435D, 0x3D7C223E7AB003EF), // 1.077596311400061984007701231389e5
    ],
    // t = 7
    [
        F64x2::from_bits(0x3EA923B08F80599A, 0x3B4E5D1474CB8ED7), // 7.492128863997167080771040272104e-7
        F64x2::from_bits(0xBEC0D878A129FEAC, 0xBB508747F0751287), // -2.008150894738791991169305312067e-6
        F64x2::from_bits(0x40F39C7CA6BEEE79, 0xBD8A40A91D3B2843), // 8.032779070943024700539121139856e4
        F64x2::from_bits(0x410994855DF32C34, 0x3DAC1E7BFE8F0561), // 2.095526708739713195059628123700e5
    ],
    // t = 7.25
    [
        F64x2::from_bits(0x3E99943724BA2234, 0x3B37DFEBB852D021), // 3.811563018337377610797492562585e-7
        F64x2::from_bits(0xBEB16EAB7B94EB67, 0xBB1B7365A2FF17A6), // -1.039046294628025735228307461362e-6
        F64x2::from_bits(0x4102F02B76056889, 0x3D9DF6CE724420CD), // 1.551414326275030975839598844092e5
        F64x2::from_bits(0x4119288C5A5C7CE4, 0xBDABAC1BB4917B74), // 4.121950882434381511883212851838e5
    ],
    // t = 7.5
    [
        F64x2::from_bits(0x3E89BBA4458FB5A6, 0xBB281E4A1994B056), // 1.917256067513430751645002898931e-7
        F64x2::from_bits(0xBEA1D396279DD5CB, 0x3B25DC3F0A523900), // -5.312713959720544684789544280410e-7
        F64x2::from_bits(0x411281F675E00F59, 0x3DB2934450297EB9), // 3.032296151125334022938305913059e5
        F64x2::from_bits(0x41290627ABB4260A, 0xBDC9E44B7CB8FDCE), // 8.199878353587996209320829744187e5
    ],
    // t = 7.75
    [
        F64x2::from_bits(0x3E7999CDA3237FD7, 0x3B177C1ED118CADC), // 9.537038961641585223672617612052e-8
        F64x2::from_bits(0xBE9204AC343B1369, 0xBB342ECD0D128AAA), // -2.684928867953261859794279547607e-7
        F64x2::from_bits(0x41224CD1420D1349, 0xBDC2AC703711BA4C), // 5.996566290060068373830522127573e5
        F64x2::from_bits(0x41392B11706CDB3C, 0xBDDDC0CE77D6C705), // 1.649425439161016497553098984985e6
    ],
    // t = 8
    [
        F64x2::from_bits(0x3E6930EBC96D9DDD, 0x3AE75760ADE60898), // 4.692207616099231625649081703488e-8
        F64x2::from_bits(0xBE8201267C1C127E, 0xBB26B6FD5D8562EB), // -1.341439297906786574291153707932e-7
        F64x2::from_bits(0x41324DE2010E4CF5, 0xBDDEECC35B3F1ADD), // 1.199586004124459930881654499598e6
        F64x2::from_bits(0x4149977328080357, 0x3DDE8A778C2504A2), // 3.354342312744538876507746496528e6
    ],
    // t = 8.25
    [
        F64x2::from_bits(0x3E58856CD2C69043, 0xBAF9B07FD03ECF9B), // 2.283713944482228170923726480449e-8
        F64x2::from_bits(0xBE71CA01D79DD13F, 0xBB0F848068A8E19A), // -6.626952666987631228217076207113e-8
        F64x2::from_bits(0x414284453A623BFD, 0xBDEB8AF0BE284988), // 2.427018456122873616239405647229e6
        F64x2::from_bits(0x415A4DD858C0D2D5, 0xBDF746A3B32C314B), // 6.895457386769015934610282971146e6
    ],
    // t = 8.5
    [
        F64x2::from_bits(0x3E479DAB884916E4, 0x3ACEFFAB07A4B0F3), // 1.099700975519550650949062908074e-8
        F64x2::from_bits(0xBE6161E55DBE1825, 0xBAFFDE70E5167398), // -3.237725440447602255894237298698e-8
        F64x2::from_bits(0x4152F0F1E2A7773F, 0x3DFBB308F003B445), // 4.965319541471301981063954390843e6
        F64x2::from_bits(0x416B5343A0FB2EFF, 0x3DFC5F7030D6A408), // 1.432630103066205833416906018644e7
    ],
    // t = 8.75
    [
        F64x2::from_bits(0x3E3681910808B5D4, 0xBAC9CCC87886536F), // 5.240114231891752419198105426726e-9
        F64x2::from_bits(0xBE50CCF2FE8C3DE6, 0xBAF20555DB6B66AB), // -1.564676202757794909372219715136e-8
        F64x2::from_bits(0x416396B9EF2E9B52, 0x3DEAF1F63CC715B0), // 1.027015947443929706747558909949e7
        F64x2::from_bits(0x417CAF66A6A03767, 0xBDDEDB1011D14475), // 3.007857041411533568005065188980e7
    ],
    // t = 9
    [
        F64x2::from_bits(0x3E253A28272EABA4, 0xBABE4FCE9760CF58), // 2.471168430872489843289241134339e-9
        F64x2::from_bits(0xBE401086AE331E68, 0xBAD771CF40379E17), // -7.480641389658946412759545273419e-9
        F64x2::from_bits(0x41747A664E4351B6, 0x3E1B6ADD2E590C7C), // 2.147286889143534909336812881111e7
        F64x2::from_bits(0x418E6D000E3F4CCE, 0x3E1EBB4CB02304E1), // 6.380748978090821385451353367766e7
    ],
    // t = 9.25
    [
        F64x2::from_bits(0x3E13D12A38B2B972, 0xBA96D34C88F70662), // 1.153504155728340160840033630029e-9
        F64x2::from_bits(0xBE2E65D2A40B2F7C, 0xBAB97ABD19EBBE10), // -3.538763310465634886516643361923e-9
        F64x2::from_bits(0x4185A2EF6A5252F2, 0xBE086A59D43A85C1), // 4.537495729019726741781245270869e7
        F64x2::from_bits(0x41A04D32C70CF077, 0xBE4FD2CC69A0E696), // 1.367473635252720863449240154964e8
    ],
    // t = 9.5
    [
        F64x2::from_bits(0x3E02508CAE8391C9, 0xBAAAB83E5774731F), // 5.330263704617491626585486669522e-10
        F64x2::from_bits(0xBE1C75FA3685ED87, 0x3AB531A57E1D7F97), // -1.656639459374066626258758935215e-9
        F64x2::from_bits(0x419719D5A65261C5, 0x3E3CA45147C26C5D), // 9.689226558045109283222473205802e7
        F64x2::from_bits(0x41B1A521CBDE3594, 0x3E1757E10A12D7EE), // 2.960347638680050386664968436801e8
    ],
    // t = 9.75
    [
        F64x2::from_bits(0x3DF0C216471C2A3A, 0x3A9B8CC305EB8444), // 2.438632135722847079048296016159e-10
        F64x2::from_bits(0xBE0A5FD1EBD1F31A, 0xBAA77CE2D76C60E3), // -7.675930651861793049430396709061e-10
        F64x2::from_bits(0x41A8EB9FC7276C0F, 0xBE3E936C524AAB5D), // 2.090475235769962895252109495492e8
        F64x2::from_bits(0x41C34A50252E268D, 0x3E6DEA9028F1BB69), // 6.472745703605515322019573414576e8
    ],
    // t = 10
    [
        F64x2::from_bits(0x3DDE5E028A1F8CDA, 0xBA4E8CCF07EBCBDA), // 1.104753255289868593355020565799e-10
        F64x2::from_bits(0xBDF831907393566E, 0x3A7D973D528B3743), // -3.520633676738923636620644825279e-10
        F64x2::from_bits(0x41BB2888418C587C, 0xBE5923E631004CA7), // 4.556411535482251409997873081065e8
        F64x2::from_bits(0x41D54C1A699EE746, 0xBE6266C862B0D319), // 1.429236134482865776118831447814e9
    ],
    // t = 10.25
    [
        F64x2::from_bits(0x3DCB3F5D60CD94D6, 0x3A6B7E2E17BF79AD), // 4.956294758320720558785060474063e-11
        F64x2::from_bits(0xBDE5F8C5EFFD582B, 0xBA8D9271FD29926C), // -1.598656693090870729415563890502e-10
        F64x2::from_bits(0x41CDE56633E7AA03, 0xBE66BB30554BF39B), // 1.003146343809875803568947583523e9
        F64x2::from_bits(0x41E7BE1BCE21E0BF, 0xBE882A9E1589DE41), // 3.186679409058684645842949506124e9
    ],
    // t = 10.5
    [
        F64x2::from_bits(0x3DB836DA085F6AFF, 0x3A5E752FCA4650FF), // 2.202274519283401643530304396356e-11
        F64x2::from_bits(0xBDD3C1E45784646B, 0xBA7D007398EE51BE), // -7.187696781451567091337852978338e-11
        F64x2::from_bits(0x41E09E71E9245FCF, 0xBE84EA3AD48D6B20), // 2.230554441136695229150651552701e9
        F64x2::from_bits(0x41FAB95DB5548865, 0xBE91D932D005E488), // 7.173692245283299180143600872290e9
    ],
    // t = 10.75
    [
        F64x2::from_bits(0x3DA550A628DC2A04, 0x3A4D5353D5788A6D), // 9.692955879668771665682310580252e-12
        F64x2::from_bits(0xBDC197E3C655ED2A, 0xBA6E0A8FC5857E59), // -3.200206141053633266681951246629e-11
        F64x2::from_bits(0x41F2A88D5E98A30D, 0x3E6681B1AF8F7703), // 5.008578025539807361562955230251e9
        F64x2::from_bits(0x420E5DEA7111EA5B, 0x3E8E99906EC1824F), // 1.630308304223943160921230999522e10
    ],
    // t = 11
    [
        F64x2::from_bits(0x3D92965D151032EA, 0x3A0BAC3C3AD0080C), // 4.226275864960359591298835450796e-12
        F64x2::from_bits(0xBDAF080902602021, 0xBA3675E7835E37EC), // -1.411144124662851733545119127150e-11
        F64x2::from_bits(0x420526DBEC13719E, 0xBEACC345C8638C56), // 1.135578253043047628513623623395e10
        F64x2::from_bits(0x42216A728E09DA9D, 0xBEC1319B9EC4AAF5), // 3.740016819692697701528300933494e10
    ],
    // t = 11.25
    [
        F64x2::from_bits(0x3D800F0851E1CB23, 0x3A27FF17E7797DB2), // 1.825665174335469456567692613172e-12
        F64x2::from_bits(0xBD9B1B5C6BF50DF1, 0xBA33624A91A35530), // -6.163390706512229683798971418205e-12
        F64x2::from_bits(0x4218356A5D4274E7, 0x3EB696F48E15D322), // 2.599384456061416379151756247713e10
        F64x2::from_bits(0x423429AA93F78782, 0x3EB6A42D77685496), // 8.659839077552932874210588710021e10
    ],
    // t = 11.5
    [
        F64x2::from_bits(0x3D6B7E7CA3A6EC4E, 0xB9EBB8B9ABB1AE84), // 7.814290183962854346130297587929e-13
        F64x2::from_bits(0xBD8774D45718AAA0, 0xBA2A0C095A1D10F7), // -2.666679967504531405901069622158e-12
        F64x2::from_bits(0x422BF863163DCAC5, 0x3EC3759BCF28A95E), // 6.006568015889603656046494363945e10
        F64x2::from_bits(0x42478EF2FABF3122, 0x3EDED5A8FFFA0450), // 2.023650727663838574492001495262e11
    ],
    // t = 11.75
    [
        F64x2::from_bits(0x3D5752B1F07BD591, 0xB9F7D144CE33F381), // 3.314401573051556788982823960803e-13
        F64x2::from_bits(0xBD741BE955EDB6A8, 0xBA0066B7B70CF2A8), // -1.143065967971401522675979369416e-12
        F64x2::from_bits(0x42404F4DD36E4E61, 0x3EEF1CED829C538F), // 1.400999543966123500409122076078e11
        F64x2::from_bits(0x425BC6F9D37570E5, 0xBEFC6F2300C3B6BB), // 4.772095134297639499338838070055e11
    ],
    // t = 12
    [
        F64x2::from_bits(0x3D439B7A11F5A8EE, 0x39E81C559B5BF38A), // 1.393184688875360839049034503196e-13
        F64x2::from_bits(0xBD6114C208E15BE4, 0x3A0BDDDB4F44AB89), // -4.854736554985308462993653997695e-13
        F64x2::from_bits(0x42533282B8F944BF, 0x3EF35DFB89997858), // 3.298072258290741761847681118236e11
        F64x2::from_bits(0x427086185756B5EF, 0xBF1C1A68C0201CA1), // 1.135507502443370742404324090458e12
    ],
    // t = 12.25
    [
        F64x2::from_bits(0x3D3056579458A40C, 0xB9C215386E78FBBD), // 5.804165630527137466681199332592e-14
        F64x2::from_bits(0xBD4CC12C98E79A59, 0xB9B01F9820B57990), // -2.043136273538622520425746726206e-13
        F64x2::from_bits(0x4266CDAB02C52D13, 0x3F0732B2663101CB), // 7.835185372574086135826771841759e11
        F64x2::from_bits(0x4283D5B9AA0304B0, 0xBF1F3B067207D917), // 2.726082986080585818364841081653e12
    ],
    // t = 12.5
    [
        F64x2::from_bits(0x3D1AFC62C7A4A98A, 0xB9BB649A86B71849), // 2.396827826078049936281668939411e-14
        F64x2::from_bits(0xBD37FC46FE0F2FBC, 0xB9CC49A5D44E82A4), // -8.521346564673856445296977249553e-14
        F64x2::from_bits(0x427B552EC71860D4, 0x3F1CBA75BEC55E99), // 1.878291935622051867402493092180e12
        F64x2::from_bits(0x429805321DEE512E, 0x3F3EB7826F3395ED), // 6.602648681364295390576690989752e12
    ],
    // t = 12.75
    [
        F64x2::from_bits(0x3D0617F77EA50C9B, 0x39AF6C2156BD8D47), // 9.811538346010540646987073539975e-15
        F64x2::from_bits(0xBD23D4511FFABB17, 0xB9C6679E314D065F), // -3.522402435158488372468500384833e-14
        F64x2::from_bits(0x4290872C7ACAA5DB, 0xBF1FB0FD580E149F), // 4.543188218537463746294533835907e12
        F64x2::from_bits(0x42AD5815B3848E19, 0x3F40A541B80F1C61), // 1.613207920698304933610638556097e13
    ],
    // t = 13
    [
        F64x2::from_bits(0x3CF1EEACDE5A021B, 0x398F45502EE38A03), // 3.981776078833335363022547078760e-15
        F64x2::from_bits(0xBD103FC396528972, 0xB9ACD1D0710F120A), // -1.443208057397262604448112020465e-14
        F64x2::from_bits(0x42A42AA66E72A6CF, 0x3F3D82032AF2872E), // 1.108670671905940474712922418564e13
        F64x2::from_bits(0x42C214632C2B5A2C, 0x3F5B18043221C9B2), // 3.975754496990834540367517231783e13
    ],
    // t = 13.25
    [
        F64x2::from_bits(0x3CDCDC68EAA89D26, 0x397AC2679B71FE6C), // 1.602105997911432307523284348109e-15
        F64x2::from_bits(0xBCFA65F79BC55CF9, 0xB99187AC4AAA6D29), // -5.861602193750886335252234564828e-15
        F64x2::from_bits(0x42B8D2A166739512, 0x3F5FE85E5E4789D0), // 2.729293004482107225999084141021e13
        F64x2::from_bits(0x42D6786CDB1FDB5B, 0xBF568922C550DC60), // 9.882591274583742049953083436462e13
    ],
    // t = 13.5
    [
        F64x2::from_bits(0x3CC707475301EA21, 0xB931142A683D3594), // 6.391673876741866650678720253362e-16
        F64x2::from_bits(0xBCE5421E416A4613, 0x3958C25AC34334F4), // -2.360142543924311285439650436375e-15
        F64x2::from_bits(0x42CED1FD356435F4, 0x3F4FCFB0768A0F5A), // 6.777449026570790722080341682030e13
        F64x2::from_bits(0x42EC2A6A485BE961, 0x3F7E218A7CEE0C01), // 2.477479786494190386062059104715e14
    ],
    // t = 13.75
    [
        F64x2::from_bits(0x3CB23871F0330D13, 0x392376201D7240FF), // 2.528600739926814900934853475644e-16
        F64x2::from_bits(0xBCD0F901792DA440, 0xB97009D41048C852), // -9.421729358681710065773523687675e-16
        F64x2::from_bits(0x42E34C6E42F8B5E5, 0x3F56A76F66CD3018), // 1.697519066003671576326930910055e14
        F64x2::from_bits(0x4301CD295D1CD00B, 0x3F74F491815DD0B2), // 6.263315179504653801160510828262e14
    ],
    // t = 14
    [
        F64x2::from_bits(0x3C9C97D1D62ADB7C, 0x39374EF3225C41EE), // 9.920205491192377266317332882133e-17
        F64x2::from_bits(0xBCBADF5CE82292BA, 0x39575E8AD7A50CF7), // -3.729310110017900679713424652832e-16
        F64x2::from_bits(0x42F85FF07445EA62, 0x3F9920670C5EFA44), // 4.288053617865341495376683521903e14
        F64x2::from_bits(0x4316B0BAC19D4E8B, 0x3FA3C2BE4A8503CC), // 1.596691411588002788595148655538e15
    ],
    // t = 14.25
    [
        F64x2::from_bits(0x3C8640181CB4ECC7, 0xB92758C8E0715E47), // 3.859823558340153006792131556253e-17
        F64x2::from_bits(0xBCA5184578AF5567, 0xB94393B93723081E), // -1.463746488588166406466940877791e-16
        F64x2::from_bits(0x430F0C04D1B0998D, 0x3FA6B889E272E8E4), // 1.092367389430577669376667874035e15
        F64x2::from_bits(0x432D298070440A58, 0xBFCA813F3A7DE1E9), // 4.104202970334507792930694999300e15
    ],
    // t = 14.5
    [
        F64x2::from_bits(0x3C712C56AE2CB3DE, 0x3912506118E1BF8D), // 1.489537454965927195298039298953e-17
        F64x2::from_bits(0xBC906BEE3445E59A, 0x39397A746F997C70), // -5.697388206185780610665177680570e-17
        F64x2::from_bits(0x4323F04FB39D3761, 0xBFAD4F86AE99389A), // 2.806124832005040442752638994888e15
        F64x2::from_bits(0x4342E4E6F4BB4C36, 0xBFEB990CD4DA6B57), // 1.063646036063652313756712367162e16
    ],
    // t = 14.75
    [
        F64x2::from_bits(0x3C5A4B5A3DE97F4E, 0x38FA9BD9F9CF3AFA), // 5.701677353184336573538811892774e-18
        F64x2::from_bits(0xBC795B474957553E, 0xB91C04CCDD1F24E4), // -2.199330691133481467552777764946e-17
        F64x2::from_bits(0x4339D2A1300CE3AC, 0x3FD5131EDFD5259B), // 7.268464155485100329292028986481e15
        F64x2::from_bits(0x4358AED28F044C39, 0xBFF0A3A8D5D6D41D), // 2.779047522737789096004406422684e16
    ],
    // t = 15
    [
        F64x2::from_bits(0x3C43F7DF8B7F2918, 0xB8EF2A3B24EFED50), // 2.164962520737992298989454038808e-18
        F64x2::from_bits(0xBC636AA0407E6A72, 0xB9033209DCA1DAC1), // -8.420567954017772766124392806839e-18
        F64x2::from_bits(0x4350DC079D73F075, 0x3FFAF4E5DCD8049F), // 1.898209956749358968478952662580e16
        F64x2::from_bits(0x437040CBB2F23E66, 0xC01C2706A0DE7717), // 7.319749203407010496188877718144e16
    ],
    // t = 15.25
    [
        F64x2::from_bits(0x3C2E162AC53F329C, 0x38C52A3FCE49B321), // 8.154986594314401873627096259974e-19
        F64x2::from_bits(0xBC4D7EB957F4A741, 0xB8ED43EBDBF63B8F), // -3.197856001111869627444267346836e-18
        F64x2::from_bits(0x436631DEBC2B8CF5, 0xC00E1C4A284ABDD1), // 4.997825757684726823618668100269e16
        F64x2::from_bits(0x43859390DA146195, 0x4000B21DA7243515), // 1.943431953365203540869706209111e17
    ],
    // t = 15.5
    [
        F64x2::from_bits(0x3C167CA306047B7B, 0xB899FDCAE49FCF37), // 3.047538152456012684173786516248e-19
        F64x2::from_bits(0xBC3638F4A15CA465, 0xB8DA5870EA6F33CF), // -1.204683204453443742272191285783e-18
        F64x2::from_bits(0x437D748F30E03A9A, 0xC0113B95B1170DE7), // 1.326549227800928276918117837600e17
        F64x2::from_bits(0x439CDF13FE614D09, 0x402B5A3ACF136106), // 5.201008840358795016762299262196e17
    ],
    // t = 15.75
    [
        F64x2::from_bits(0x3C00ACCAC8666808, 0x389766ED132066B6), // 1.129939865914347374486471568287e-19
        F64x2::from_bits(0xBC209C15C761B3B0, 0x38C7514A87059795), // -4.502070189805177478784587153769e-19
        F64x2::from_bits(0x4393B3D96DB1E68A, 0x401BB66B7E1F49CB), // 3.549293431777941829281444269749e17
        F64x2::from_bits(0x43B37811978CD8C2, 0xC045291B01664A60), // 1.402890626526659029678863358549e18
    ],
    // t = 16
    [
        F64x2::from_bits(0x3BE889B6799D2C7B, 0xB847EC36F24EC651), // 4.156888828917024394747937619181e-20
        F64x2::from_bits(0xBC08A2043B3C677A, 0xB87E14DA835D1B23), // -1.669188676838180955915934128152e-19
        F64x2::from_bits(0x43AA91673209379F, 0xC00B53302731E936), // 9.572123906049186525843808114904e17
        F64x2::from_bits(0x43CA7691C1BF39BC, 0xBFBF7371B233F4FA), // 3.813743507121862655877144712442e18
    ],
    // t = 16.25
    [
        F64x2::from_bits(0x3BD1EA36171C4159, 0xB861F6A2E624FB76), // 1.517450025709020168768716499909e-20
        F64x2::from_bits(0xBBF21F607BC4D99B, 0x389E0E57FAB2580A), // -6.140164186475070829665663562191e-20
        F64x2::from_bits(0x43C20DF2C95EAD6E, 0x406E0CDAA1B8D76B), // 2.601925628141427952401688443202e18
        F64x2::from_bits(0x43E21FF78FA7E331, 0xC0877CA7073E5D0E), // 1.044827690663841102441844321517e19
    ],
    // t = 16.5
    [
        F64x2::from_bits(0x3BB9F55B9BBD01DD, 0x384B10BCB097597E), // 5.496911172967060763600712738318e-21
        F64x2::from_bits(0xBBDA7555BF8C5A2B, 0x384BCF6FA91F213A), // -2.241108542525297257484245606679e-20
        F64x2::from_bits(0x43D8BB092FE10979, 0xC079E6B0975454BF), // 7.128112715195081313581886931021e18
        F64x2::from_bits(0x43F904FCC1717DDD, 0xC09424ADFB0B7C22), // 2.884549853788804171883009702735e19
    ],
    // t = 16.75
    [
        F64x2::from_bits(0x3BA2A9E31AF21570, 0x382FB0024EE0C69C), // 1.976087801413193513358794026349e-21
        F64x2::from_bits(0xBBC32A415B2E1578, 0xB8659D72B1ADC19F), // -8.116718524732295808677260375320e-21
        F64x2::from_bits(0x43F111CBC800798A, 0x4080F6564D734A70), // 1.967981172977150415879213991233e19
        F64x2::from_bits(0x44116695405CD7B1, 0x403BC869DC25DBE2), // 8.024663359599019625178286529469e19
    ],
    // t = 17
    [
        F64x2::from_bits(0x3B8AA2884DD9FB25, 0x3822360C558684C9), // 7.050197298388614542441447162870e-22
        F64x2::from_bits(0xBBAB8D3B4A648A48, 0x384FD6D0C4D1CEB9), // -2.917148219293313793259713373315e-21
        F64x2::from_bits(0x4407BECEF08F682A, 0xC07FFEBD01BACE03), // 5.475303811330586982407885577230e19
        F64x2::from_bits(0x442863570960CE96, 0x40A9E6C655455069), // 2.249400291065726927233873693143e20
    ],
    // t = 17.25
    [
        F64x2::from_bits(0x3B72DCE3AA2A38A3, 0x3816A50FC1DF6CD9), // 2.496477263592584405057664688714e-22
        F64x2::from_bits(0xBB93A7513DA528B1, 0x383A75B17AA317D0), // -1.040452041235295247058721086534e-21
        F64x2::from_bits(0x4420A4822C16EE77, 0x40A23B64ABCAFDFA), // 1.535009954127735380776966231761e20
        F64x2::from_bits(0x4441383B217598A3, 0xC0EA44B890466298), // 6.352930366635359062662323883128e20
    ],
    // t = 17.5
    [
        F64x2::from_bits(0x3B5A84B97BBA35C1, 0xB7BF13C0C24F450A), // 8.774220823294709737521657693715e-23
        F64x2::from_bits(0xBB7BD3DC52C12E20, 0xB80979FABE453329), // -3.682949628790096691918044688854e-22
        F64x2::from_bits(0x443781A13458FEE1, 0xC0DED0A409027412), // 4.336159183543006998374369500987e20
        F64x2::from_bits(0x44587FB579A80632, 0xC0F3821D326FDD3C), // 1.807697011780784408542175186287e21
    ],
    // t = 17.75
    [
        F64x2::from_bits(0x3B4280B50513C75D, 0x37AF118BEED8FD54), // 3.061025224748163232666220888200e-23
        F64x2::from_bits(0xBB638D913B1490B3, 0x37EB709DF0B7EFBE), // -1.293908822436053441068011933002e-22
        F64x2::from_bits(0x4450B9CCF38B674D, 0x40C7EB0B54536EFE), // 1.234144995006710296534088510923e21
        F64x2::from_bits(0x44718EAE66EBBB37, 0x41157E0E802A0856), // 5.182014678640641793923625160341e21
    ],
    // t = 18
    [
        F64x2::from_bits(0x3B29A15F044779A9, 0x3791DC1C171677C0), // 1.060046682524795565634125585294e-23
        F64x2::from_bits(0xBB4B45FC155A9C58, 0xB7C9E8A59C27A1BF), // -4.512001860681941889168941687420e-23
        F64x2::from_bits(0x4467FB0B68D69E46, 0xC0F40DB726F0D698), // 3.538918250356568669988552992975e21
        F64x2::from_bits(0x448959F229E6CD81, 0xC125612FE4E7B7D1), // 1.496479650328785068426405291963e22
    ],
    // t = 18.25
    [
        F64x2::from_bits(0x3B119F53AD46603F, 0xB7BF3006E2A2142E), // 3.644220818093478600760544747445e-24
        F64x2::from_bits(0xBB32E16D6AEF4E3E, 0xB7D74434626ED5AC), // -1.561764628191697478015882815802e-23
        F64x2::from_bits(0x448151B0F62B5B60, 0xC111EB89FE22F7BD), // 1.022339598822760484175750181973e22
        F64x2::from_bits(0x44A26FE05700250E, 0xC135FAD041D86A72), // 4.353317533669099477534374279151e22
    ],
    // t = 18.5
    [
        F64x2::from_bits(0x3AF80EABDB31A22D, 0x3779C487C73E19A6), // 1.243733766971940457467866405745e-24
        F64x2::from_bits(0xBB19F301F39C9C11, 0x379A1CA56174AE2F), // -5.366178823414727709383683897215e-24
        F64x2::from_bits(0x4499337518772B08, 0xC132A7E927C62206), // 2.975209591110722318133484463298e22
        F64x2::from_bits(0x44BB0338CE46E91D, 0x4155D89B377D9F42), // 1.275633285567799093745728670424e23
    ],
    // t = 18.75
    [
        F64x2::from_bits(0x3AE04D885F453106, 0x37865F7833AB6615), // 4.214191086870341438613927676805e-25
        F64x2::from_bits(0xBB01B3C9CD456C73, 0xB794EEB10299EE3D), // -1.830374713082134616442949605010e-24
        F64x2::from_bits(0x44B278333A7112C1, 0xC14B6E988A9E6572), // 8.721989736766870134036691704113e22
        F64x2::from_bits(0x44D3EE8ADECABA50, 0x416A99B96AB807A8), // 3.765011796895590182291953349646e23
    ],
    // t = 19
    [
        F64x2::from_bits(0x3AC5F01B0B5B7CEC, 0xB7639392ABDD68F1), // 1.417704377793352718937684362985e-25
        F64x2::from_bits(0xBAE7FA58A7D63BBE, 0xB76B1EB8D5AA3386), // -6.198145827130015058571783813917e-25
        F64x2::from_bits(0x44CB4500CB5540AD, 0x4168AADD5BBEEA15), // 2.575535552234458545743148670550e23
        F64x2::from_bits(0x44EDA03B3D9A14AA, 0xC16052C91009C06A), // 1.119235006339588781074871498810e24
    ],
    // t = 19.25
    [
        F64x2::from_bits(0x3AAD4FA200AEABF9, 0xB73DFFCAA9B8F9AB), // 4.735458064902290919627065233421e-26
        F64x2::from_bits(0xBAD01F597D0AB6A5, 0x377F77FA9DDF02DF), // -2.083779016580934019645703414195e-25
        F64x2::from_bits(0x44E446E6240A5152, 0x417CFA043BB2CD0C), // 7.660414810689987684352677311526e23
        F64x2::from_bits(0x45062CC37EFB1AA8, 0xC191169DF6CCECFB), // 3.350969802724856976529538299877e24
    ],
    // t = 19.5
    [
        F64x2::from_bits(0x3A937166718B66C3, 0x3705D7DB971349F1), // 1.570590561517818377578549910856e-26
        F64x2::from_bits(0xBAB586BD89CEEB4F, 0x375928B29FDDF0E4), // -6.955532236463624257855370790846e-26
        F64x2::from_bits(0x44FE5F28C56131E1, 0xC19C4546AC12DBE4), // 2.294821968548648666967636981583e24
        F64x2::from_bits(0x4520B7382E5DD354, 0x41C2DFD0E63A4A4D), // 1.010402126946309859020846045539e25
    ],
    // t = 19.75
    [
        F64x2::from_bits(0x3A799D098B2FA5B1, 0x371E223CC9701AED), // 5.172590364317375193053113209864e-27
        F64x2::from_bits(0xBA9C89A1F1686CA2, 0xB732FB77EFD3C771), // -2.305251175214061412980809946465e-26
        F64x2::from_bits(0x4516E896E31C86D4, 0x41BB91E9959A6C70), // 6.923685108065512216193429603217e24
        F64x2::from_bits(0x453961020F6D20CA, 0x41D7DA5E030889D8), // 3.068125304427296895909684413341e25
    ],
    // t = 20
    [
        F64x2::from_bits(0x3A60C0E67717E4A7, 0xB7057DCF2EAABB4A), // 1.691672868670540313553560212509e-27
        F64x2::from_bits(0xBA82C8726B3F76B1, 0xB725FFBCD3C4D8E4), // -7.586391625748354960515371705913e-27
        F64x2::from_bits(0x453166E55190C767, 0xC1823ACA6ED91722), // 2.103765049651103814494789014400e25
        F64x2::from_bits(0x455366B452E27AB8, 0xC1F69AEEE668B9B1), // 9.381839336133964349106216945467e25
    ],
    // t = 20.25
    [
        F64x2::from_bits(0x3A45C3CEE5D250E6, 0x36ED9ED4F2FB9B18), // 5.494227069840214100466819558336e-28
        F64x2::from_bits(0xBA688D5D775F0E27, 0xB6B1269A90172D04), // -2.479139402284829034933756199359e-27
        F64x2::from_bits(0x454A9FD78441B0EB, 0x41E89D98BA91CD9A), // 6.437380632165058238668949255635e25
        F64x2::from_bits(0x456DDEA38A0679AB, 0x4207B0AE4B7F43F5), // 2.888818471446074079307226719082e26
    ],
    // t = 20.5
    [
        F64x2::from_bits(0x3A2C14A39AE69616, 0xB6C631D5BD272806), // 1.772136354314942104047501740834e-28
        F64x2::from_bits(0xBA4FDEC21E097596, 0x36D1CBAEC9FCE314), // -8.045156793755489589419322092101e-28
        F64x2::from_bits(0x4564828B2B06B8C2, 0xC1EF95345E9A703B), // 1.983599298628833545923043951738e26
        F64x2::from_bits(0x45872718319CC855, 0x422B0A0DAB75C9B2), // 8.956771981049161487518816588939e26
    ],
    // t = 20.75
    [
        F64x2::from_bits(0x3A11FD97DAF0516E, 0x36B993050DC39AE5), // 5.676831888432148240329368359124e-29
        F64x2::from_bits(0xBA348AA69BA20509, 0xB6DCE553A7E5E0D5), // -2.592715025005044673274192339155e-28
        F64x2::from_bits(0x457FD1C91DC36549, 0xC21016B7C356D3EC), // 6.154781590912205851398018982930e26
        F64x2::from_bits(0x45A211E09276C6BB, 0xC24B1EED871B0B13), // 2.796171213796595389872533961913e27
    ],
    // t = 21
    [
        F64x2::from_bits(0x39F6E5411D1E83DD, 0xB6843E4C5C266643), // 1.806138442470383379480240357006e-29
        F64x2::from_bits(0xBA1A4C3E3A148266, 0x36A8DCD6470AE1F3), // -8.298130258300445749861472396706e-29
        F64x2::from_bits(0x4598DA7D73F82633, 0xC1FDF57C17079D0E), // 1.922944624863808351294209679524e27
        F64x2::from_bits(0x45BC66164EB91385, 0xC24A4599D6592750), // 8.788996052843163912600572749693e27
    ],
    // t = 21.25
    [
        F64x2::from_bits(0x39DCF0D2224FFDD9, 0x36702083323A7B75), // 5.707547960312372586262917125723e-30
        F64x2::from_bits(0xBA00B7F2EC328AA8, 0xB6948D72C415140F), // -2.637722169248672009473273123416e-29
        F64x2::from_bits(0x45B38BC89794A303, 0xC22FE0182600C3A6), // 6.049203144997326389140253951618e27
        F64x2::from_bits(0x45D677C363FBA3ED, 0xC268A1AC169B607C), // 2.781382037771260863677512374898e28
    ],
    // t = 21.5
    [
        F64x2::from_bits(0x39C22B054E3CB54E, 0x365DD8329D47C7C6), // 1.791508017269441364627314273605e-30
        F64x2::from_bits(0xB9E51CE92EB7DD28, 0xB6229FC5C83EA049), // -8.327583751655689469681790273131e-30
        F64x2::from_bits(0x45CEF445698ED73D, 0xC26BBAAA35FD185D), // 1.915971197022905878519783835924e28
        F64x2::from_bits(0x45F1E55357AA44E6, 0xC238761C8BD24FC9), // 8.861572405133780813735811178969e28
    ],
    // t = 21.75
    [
        F64x2::from_bits(0x39A6A87E4E25DF63, 0xB63D339B2BFA3EC2), // 5.585671946326008255313642866168e-31
        F64x2::from_bits(0xB9CA7B7DEEC1B62A, 0xB630068D762636C0), // -2.611365300982023152831268980779e-30
        F64x2::from_bits(0x45E8AD4EC76D7E3A, 0xC28360059D5D35B9), // 6.109725140706523226527212040928e28
        F64x2::from_bits(0x460CB3394777D8E9, 0x42A1CCE5EA5815FB), // 2.842319516784941480711063582830e29
    ],
    // t = 22
    [
        F64x2::from_bits(0x398C11F944369D4C, 0xB62106E99F60E9F9), // 1.729960240353698285290118861065e-31
        F64x2::from_bits(0xB9B07F4C042E2834, 0x365296F8290A2CB2), // -8.133774038447540275371838871712e-31
        F64x2::from_bits(0x4603CE401E00EF1A, 0xC2A8C14DD6FD24EF), // 1.961458140918622869051977567815e29
        F64x2::from_bits(0x46272AE72B233863, 0x42BB2E233DBAA953), // 9.177627869081825969553264954013e29
    ],
    // t = 22.25
    [
        F64x2::from_bits(0x397145C651BC32CD, 0x3618F7ECC71BF969), // 5.322517670348446414297405391710e-32
        F64x2::from_bits(0xB9946AB7143832EC, 0x36313F5D7DF6095D), // -2.516571805637467864957473375462e-31
        F64x2::from_bits(0x462000B424360C86, 0xC2A39F937724C986), // 6.339341889389899926174666868752e29
        F64x2::from_bits(0x4642D373766500F0, 0xC2EB7C6C9795178E), // 2.983095695801502263424349067076e30
    ],
    // t = 22.5
    [
        F64x2::from_bits(0x39551DF91D7D7F8E, 0xB5E710A2C1A55067), // 1.626800856851621560690369682322e-32
        F64x2::from_bits(0xB97919A0177A9F53, 0x35E16AE4FC96F188), // -7.734565180555525410121522783927e-32
        F64x2::from_bits(0x463A08657A47E658, 0xC2CCC8B93FBF2D9A), // 2.062530784412252307988500480421e30
        F64x2::from_bits(0x465ECC61382C9EB7, 0xC2ED508A1FC91CAB), // 9.760389394046413864458797561627e30
    ],
    // t = 22.75
    [
        F64x2::from_bits(0x3939A608E6033F4A, 0x35D00F33D0FB2AEE), // 4.939735531102591128012267120839e-33
        F64x2::from_bits(0xB95EA76EFA71E131, 0x35F9B6ACEFD72E89), // -2.361501122115945893149010060122e-32
        F64x2::from_bits(0x465550B9063F61F8, 0x42E198B3DE301FFC), // 6.755095577461162425329212555520e30
        F64x2::from_bits(0x46795BA60BDD0371, 0xC302B09E14E47FD8), // 3.214508698335066495594230365799e31
    ],
    // t = 23
    [
        F64x2::from_bits(0x391EF330351063F6, 0xB5B96A9F28B11989), // 1.490186095767206091032867927397e-33
        F64x2::from_bits(0xB942987D531D3818, 0xB5DC4AEDC7C800A5), // -7.162788572866303810261025304351e-33
        F64x2::from_bits(0x46719166D91641D4, 0x42F6065E378AFFCA), // 2.227005479622860894543865975193e31
        F64x2::from_bits(0x469503E7637E6761, 0xC305781DD2CE4684), // 1.065599743705194806302759226059e32
    ],
    // t = 23.25
    [
        F64x2::from_bits(0x39028D848C698F5E, 0x35AC13D048C4A703), // 4.466425043882071934639960869128e-34
        F64x2::from_bits(0xB9266A089A4AD0D8, 0x35BE7ECE8D1D0329), // -2.158408347743438251061010402359e-33
        F64x2::from_bits(0x468D26284D3B0A6C, 0x432019BD624A468D), // 7.390162768099166769125324050516e31
        F64x2::from_bits(0x46B18790E8906B23, 0x43455C88CA7F6258), // 3.555416121957056162564379528900e32
    ],
    // t = 23.5
    [
        F64x2::from_bits(0x38E619875F7990CB, 0x3558E0A37D759685), // 1.330078628939432760230020425566e-34
        F64x2::from_bits(0xB90AD76FEEC3E7ED, 0x35AC2ADB90BD6F15), // -6.461868523316313294311633395756e-34
        F64x2::from_bits(0x46A85718E33C7A8A, 0x4349A01E74EF5B4A), // 2.468391914891901584483506380284e32
        F64x2::from_bits(0x46CD6EED49EF6C80, 0xC36009B26283C283), // 1.193956827424054106319549569822e33
    ],
    // t = 23.75
    [
        F64x2::from_bits(0x38CA28042F560086, 0xB564A9C2D2E51EED), // 3.935564694626075600511541604654e-35
        F64x2::from_bits(0xB8EFEF9D2AEDEDAC, 0x352EB159F628E7A9), // -1.922077599325106100464611595483e-34
        F64x2::from_bits(0x46C474EE0A572F88, 0xC36143FA1955D615), // 8.298246577160128185741856190831e32
        F64x2::from_bits(0x46E8DE8D90013BF6, 0xC386BBFE9EC14082), // 4.035282352346704257892295919530e33
    ],
    // t = 24
    [
        F64x2::from_bits(0x38AEC2AF9130B3C3, 0x3535AF4DACD2E25D), // 1.157081085398542465860530587135e-35
        F64x2::from_bits(0xB8D2E070B1305EF0, 0x3573D0153CC979A8), // -5.680506160122678296411564484037e-35
        F64x2::from_bits(0x46E14DD1DB950C31, 0x438708D7B1EE6DE2), // 2.807731836810731452464929983807e33
        F64x2::from_bits(0x470525CD09088CE1, 0xC3852B2D83593058), // 1.372561596919796143033889920817e34
    ],
    // t = 24.25
    [
        F64x2::from_bits(0x3891F9187FA40A16, 0xB531249BDDA9E1D1), // 3.380351250108078482829158206209e-36
        F64x2::from_bits(0xB8B62C3D462C53CF, 0xB5569D43AC3B8B4C), // -1.668096823049951117371409043549e-35
        F64x2::from_bits(0x46FD765DC9758AF6, 0xC38D57DD3D94D23D), // 9.561085237402640009434457889421e33
        F64x2::from_bits(0x472218F5B4F0261C, 0xC3C17182D6C0BD98), // 4.698379409077722906176102473746e34
    ],
    // t = 24.5
    [
        F64x2::from_bits(0x3874DEDEEC119C60, 0xB4F434A4098E2574), // 9.813303797462994780418646960982e-37
        F64x2::from_bits(0xB899E106925A4394, 0x352B7CC2958E642A), // -4.867300156198381662215621316097e-36
        F64x2::from_bits(0x47193DFFCEE5F205, 0x43BE6A0F5E9FECD3), // 3.276622891563338239823571907070e34
        F64x2::from_bits(0x473F2BBED18F8E4F, 0x43D661D03B385843), // 1.618484644343298731873749168680e35
    ],
    // t = 24.75
    [
        F64x2::from_bits(0x3858155819F16C61, 0x34D36E9099CE5E13), // 2.830987193098951230827192244980e-37
        F64x2::from_bits(0xB87E038A7A7D8F22, 0xB5168A80F66053BF), // -1.411243587399390142045911061222e-36
        F64x2::from_bits(0x4735C3985D06AED1, 0x43AB923D94FC3542), // 1.130053753728306405447932934710e35
        F64x2::from_bits(0x475B037127866873, 0xC3F2E838FF4B7C45), // 5.610473096960803910109679873265e35
    ],
    // t = 25
    [
        F64x2::from_bits(0x383B9E0E76611475, 0xB4D129DA5BBE1B68), // 8.116026824691386683758343296410e-38
        F64x2::from_bits(0xB8614B94473EEF4D, 0x35058225AA92B442), // -4.066089337243281005322614298216e-37
        F64x2::from_bits(0x4752E246B32FCD6D, 0x43E94B7436FD396B), // 3.922030778041381773803850112158e35
        F64x2::from_bits(0x47778EB2C7842322, 0xC406E2B6EB3DE343), // 1.957073508323330897013266831871e36
    ],
    // t = 25.25
    [
        F64x2::from_bits(0x381F79365EF75549, 0xB4A5EAE2ED41DD4D), // 2.312306553658810288991498685193e-38
        F64x2::from_bits(0xB843CECC119C362B, 0xB4D7B39394147F4C), // -1.164197929022911821060748735476e-37
        F64x2::from_bits(0x47707CEF7C4ACD58, 0xC3CF60590F1BA39C), // 1.369771880397158320381871762339e36
        F64x2::from_bits(0x4794ABFC115E3E35, 0xC41223D8FAA6CE86), // 6.869388805443362555204270532057e36
    ],
    // t = 25.5
    [
        F64x2::from_bits(0x3801D2BF30AE904C, 0x34A2F90E759A4343), // 6.547220618442566694960128709936e-39
        F64x2::from_bits(0xB8268B50A87B540B, 0xB4BFDE55B2B29F22), // -3.312572393834598873189071172136e-38
        F64x2::from_bits(0x478CF8FD2BD758E4, 0xC4093E37AB818642), // 4.813901053247716562098433125896e36
        F64x2::from_bits(0x47B2409C5689A2A8, 0x4453A287C51D9D9F), // 2.426158183519094770248334805797e37
    ],
    // t = 25.75
    [
        F64x2::from_bits(0x37E40FF2B77C3894, 0xB47252EA7FAFA1C4), // 1.842431027856755412878669954997e-39
        F64x2::from_bits(0xB8097FEDF3AAF76A, 0xB4906826D789888E), // -9.367119444609421763379351184299e-39
        F64x2::from_bits(0x47A99D24E87BA6C0, 0x444272EB0FC7F8F6), // 1.702331954368027273973536158733e37
        F64x2::from_bits(0x47D0373CE4D93B50, 0x446C9DB4D3F8240C), // 8.621783734242117775271128321829e37
    ],
    // t = 26
    [
        F64x2::from_bits(0x37C671CCA74B6A23, 0x34483AB6F91D5220), // 5.153011145198340303609474433999e-40
        F64x2::from_bits(0xB7ECAA3BC53E4097, 0x34714A76815571A8), // -2.632462157874610049373415701569e-39
        F64x2::from_bits(0x47C6C8EDF3D7D9B6, 0xC46C3576DDAEE1FF), // 6.057260305859428509082886910131e37
        F64x2::from_bits(0x47ECFD778DF5637B, 0xC487B19048418587), // 3.082756787283583604523677442580e38
    ],
    // t = 26.25
    [
        F64x2::from_bits(0x37A8F4EFBFF4C5C2, 0xB446AF7DB4854C0C), // 1.432449110083830188271857053669e-40
        F64x2::from_bits(0xB7D003440B0BF166, 0xB47358728F00F2D5), // -7.352697422557221261871394161113e-40
        F64x2::from_bits(0x47E464B8B4960D1C, 0xC467040BB409EB16), // 2.168602869794137920513567724626e38
        F64x2::from_bits(0x480A1292212A2AF0, 0xC488855E81273A86), // 1.109003298983782845438772703492e39
    ],
    // t = 26.5
    [
        F64x2::from_bits(0x378B950548390164, 0x342C38CCB6F4758F), // 3.957838947992781643763624430183e-41
        F64x2::from_bits(0xB7B1C7E47A633CD4, 0xB458BD59F3187F13), // -2.041139187038112398802325818618e-40
        F64x2::from_bits(0x48025D771B1D1656, 0xC49FC1F158B61AF5), // 7.811649051800804390514961193881e38
        F64x2::from_bits(0x4827976DBEF4D233, 0xC49309EEE7DB2197), // 4.013888849978618531150331487172e39
    ],
    // t = 26.75
    [
        F64x2::from_bits(0x376E4CBEC7707A07, 0xB40CB86BA81EB9EB), // 1.086951517994317993599410320333e-41
        F64x2::from_bits(0xB7939FCBC066715A, 0x3439845ADE15C06B), // -5.631869988894562477434379013447e-41
        F64x2::from_bits(0x4820A3BB043213EC, 0x44C6FAFD00A74B54), // 2.831076539328857605775850121989e39
        F64x2::from_bits(0x484579E055C593CC, 0x44D57821B119B6C1), // 1.461586221396899960105096850109e40
    ],
    // t = 27
    [
        F64x2::from_bits(0x37508AEFB6988C77, 0x33E912FED30C5ACE), // 2.967204922863166517969256954606e-42
        F64x2::from_bits(0xB775871A86746503, 0x340B09DBF6CD600C), // -1.544540206238969500355727403726e-41
        F64x2::from_bits(0x483E55F39952F8CD, 0xC4CEBB4E7DCCCF95), // 1.032272022518887422306337679677e40
        F64x2::from_bits(0x4863AB1ACA90ED00, 0xC50E6B437107F657), // 5.354241653619452830407629014963e40
    ],
    // t = 27.25
    [
        F64x2::from_bits(0x3731F4A6C08E67A7, 0xB3CC4C4D7D6D7121), // 8.051600642213282321890311170986e-43
        F64x2::from_bits(0xB757794D97C214B9, 0xB3DCBDC4CF361690), // -4.210413600670857774764194386969e-42
        F64x2::from_bits(0x485BD1F195B9B498, 0x44BEE79F218C9F2D), // 3.786674775327769971449047930715e40
        F64x2::from_bits(0x48821F004530E804, 0xC4E01C99BBFF8BBD), // 1.973212824258318248323576432572e41
    ],
    // t = 27.5
    [
        F64x2::from_bits(0x37135F95F32953AB, 0x33AE177673BC10B6), // 2.171831213559361362676085233251e-43
        F64x2::from_bits(0xB739714B91FAFE2E, 0xB3BBA26E9EEE5747), // -1.140883819522246484852696456633e-42
        F64x2::from_bits(0x4879AAB2B8A358BC, 0x45009F64EE3982B7), // 1.397432944968437785910886313503e41
        F64x2::from_bits(0x48A0CBA12A9801F1, 0x45217331073304D9), // 7.315440611299919896239979748539e41
    ],
    // t = 27.75
    [
        F64x2::from_bits(0x36F4C77DB44046DB, 0x338B8B738D4BEEFE), // 5.823589541921340238182606901258e-44
        F64x2::from_bits(0xB71B697836C3BBBF, 0x33A1F8E57436593D), // -3.072990423064113499656218628441e-43
        F64x2::from_bits(0x4897D27AE0538491, 0xC53A9D54806928C9), // 5.188013010169519207296926337843e41
        F64x2::from_bits(0x48BF51A3B8EC4AE4, 0xC55CD8F9F12B9512), // 2.728261359943680546853271337927e42
    ],
    // t = 28
    [
        F64x2::from_bits(0x36D627E1983308E3, 0x336D98F8B6D3E035), // 1.552343448341592586671362486312e-44
        F64x2::from_bits(0xB6FD5BCD19B9E350, 0x339ED94D65884984), // -8.228031752356181205744313508361e-44
        F64x2::from_bits(0x48B63DFE470F0E81, 0xC54F751B9E8CD918), // 1.937565507750294182938258996598e42
        F64x2::from_bits(0x48DD5FAE5D34D49B, 0x4552B1BE630B69E1), // 1.023525954601976155400574895186e43
    ],
    // t = 28.25
    [
        F64x2::from_bits(0x36B77C1FB39CC4F5, 0xB35663965AB6E491), // 4.113662058837776934271624808707e-45
        F64x2::from_bits(0xB6DF41F7289BEC5A, 0x3370A93443C51998), // -2.190066780358608210889221669946e-44
        F64x2::from_bits(0x48D4E3EE9114E456, 0xC5759C5D2AE8AC61), // 7.279227701984151944000605970754e42
        F64x2::from_bits(0x48FBB647D0EC715B, 0xC577D948027B0AF3), // 3.862493713469142955537162368200e43
    ],
    // t = 28.5
    [
        F64x2::from_bits(0x3698BF8A603F3CD5, 0x3324BA12C4B62CD8), // 1.083738185658727013848100135431e-45
        F64x2::from_bits(0xB6C08ABC27BE8A61, 0xB33E23FD2E9AFBF2), // -5.795046510086222979055271995667e-45
        F64x2::from_bits(0x48F3BC9F4B6890D1, 0x455ED602428BA3DC), // 2.750909210366431128034101580433e43
        F64x2::from_bits(0x491A4C4773321970, 0x45B1171068BBBFC1), // 1.466160555155029200497123621126e44
    ],
    // t = 28.75
    [
        F64x2::from_bits(0x3679ED83AE824BA7, 0xB3093051AA26DA54), // 2.838482272275884300484188252299e-46
        F64x2::from_bits(0xB6A167E6224B79B0, 0x3325831F6EA52B79), // -1.524424895251976954451674772871e-45
        F64x2::from_bits(0x4912C1BCEC4B8F37, 0xC5B6B891680B2A36), // 1.045725929833740160517235224117e44
        F64x2::from_bits(0x49391A21B0548642, 0x45DD7D5DC14E91C6), // 5.597950130986030386347289196625e44
    ],
    // t = 29
    [
        F64x2::from_bits(0x365B0199A20F5748, 0x3266D0F2B14F4CDB), // 7.391370577886336605077244315661e-47
        F64x2::from_bits(0xB6823547ACC9B91F, 0x33045EC169CE6440), // -3.986721483253100799124307423348e-46
        F64x2::from_bits(0x4931EE12839A3BCA, 0x45C22C9833B93093), // 3.998516924390779505434155079650e44
        F64x2::from_bits(0x495819A44065D6D1, 0xC5EA9060FAE56E7A), // 2.149806335282298745266016507796e45
    ],
    // t = 29.25
    [
        F64x2::from_bits(0x363BF7A2468CC099, 0x32C249CF17B49E95), // 1.913601679571197951653452751060e-47
        F64x2::from_bits(0xB662EFD306EDB49A, 0xB2FCA77F0CE2CCE1), // -1.036567502013828289195908394065e-46
        F64x2::from_bits(0x49513D5B246BC29B, 0xC5E00CA8047E8D42), // 1.537830127695913775464404546057e45
        F64x2::from_bits(0x497745C00B9CF09C, 0x460059A3885BE33E), // 8.303891791141296336686616124751e45
    ],
    // t = 29.5
    [
        F64x2::from_bits(0x361CCBD6A28A7A21, 0xB2BBE98C13389621), // 4.925797392817714742351921615737e-48
        F64x2::from_bits(0xB64394B20468F869, 0x32E8A4052D287A77), // -2.679551062709984643564851811720e-47
        F64x2::from_bits(0x4970AC1CBA86D587, 0xC61625C7F30E50C1), // 5.948880195192517370126807163517e45
        F64x2::from_bits(0x49969A5E17EF9FD8, 0x46250B7EE4A498D3), // 3.226007710464641516749243320021e46
    ],
    // t = 29.75
    [
        F64x2::from_bits(0x35FD7AEB92FE0911, 0xB29222CE7AA0824B), // 1.260696375727886170752018033456e-48
        F64x2::from_bits(0xB62421588EBBE178, 0xB2C1219CBEC07EA2), // -6.886840567556511320980774133732e-48
        F64x2::from_bits(0x4990378A819E0304, 0x463FB5214853D770), // 2.314561473025455594459118943679e46
        F64x2::from_bits(0x49B6143D78A1C843, 0x465F4C1E8E2B9FFB), // 1.260491667604664854166661564409e47
    ],
    // t = 30
    [
        F64x2::from_bits(0x35DE0227A9491606, 0x327CF48E6CA75814), // 3.208217591550495571075286933185e-49
        F64x2::from_bits(0xB604939581BB059E, 0xB2A788733A553D96), // -1.759876581432725982082104692405e-48
        F64x2::from_bits(0x49AFBADB75FD756F, 0xC602EE62D7039340), // 9.057288512151306951894813687353e46
        F64x2::from_bits(0x49D5B0D8AE0CE676, 0xC671CA07D4EEDC1C), // 4.953304512891299042078921677826e47
    ],
    // t = 30.25
    [
        F64x2::from_bits(0x35BE5F75401E262C, 0x3248786450A9B32F), // 8.117956942534393861793209113107e-50
        F64x2::from_bits(0xB5E4E9A145ECF087, 0x325490202EF8A262), // -4.471560394793337836456022371692e-49
        F64x2::from_bits(0x49CF38271EBDC62F, 0xC647C84038C6A3FD), // 3.564619269330644075905891419690e47
        F64x2::from_bits(0x49F56E514159BE9B, 0xC68FBBCA1DE2BAC8), // 1.957583463871108604296463843119e48
    ],
    // t = 30.5
    [
        F64x2::from_bits(0x359E91701F79ABE7, 0xB209D0B372D8D697), // 2.042534616669260159633022584500e-50
        F64x2::from_bits(0xB5C52229A0640D4A, 0xB23E0418952F8607), // -1.129694662719963766022217213972e-49
        F64x2::from_bits(0x49EEE47E9617E845, 0x468C76E3C5ECBD08), // 1.410922562017126984130718566720e48
        F64x2::from_bits(0x4A154B60B33E6A77, 0xC6B31108EBADB945), // 7.780466294409502806154116948406e48
    ],
    // t = 30.75
    [
        F64x2::from_bits(0x357E976E33950CCB, 0x321C6B2ADCAA97B9), // 5.110246826661177406994226075271e-51
        F64x2::from_bits(0xB5A53C5A4A9D68A4, 0x32473A126F65C5A1), // -2.837908525573153107938792091190e-50
        F64x2::from_bits(0x4A0EBE40AE0562AA, 0x46619398E4745D0D), // 5.616399971680829262457510921850e48
        F64x2::from_bits(0x4A35474E12AE41E3, 0x46D0182DAD31483B), // 3.109861381412079766988027675228e49
    ],
    // t = 31
    [
        F64x2::from_bits(0x355E71830C073025, 0xB1EAD252BD990913), // 1.271375950934322052293414993223e-51
        F64x2::from_bits(0xB58537E204E4D9A7, 0xB221E26BEFE9B2D2), // -7.088937934666453203015320959552e-51
        F64x2::from_bits(0x4A2EC496F9335C9F, 0x46ABE247F664F1C9), // 2.248368906561252626198318140906e49
        F64x2::from_bits(0x4A5561E7C471C935, 0x46F5B413B5E415CF), // 1.250019002943485446273017395541e50
    ],
];

// GENERATE: bessel_table j0 15 0 1 11 10
// J0_TBL[i] has the coefficients `K0, ..., K20` such that
// f(c + d) ~= K0 + K1 * d + ... + K20 * d^20, where `c` is very close to the
// (i + 1)-th positive zero of j0, for the values of the interval of that
// zero, [(i + 0.25) * π, (i + 1.25) * π) (the first one starts at 1)
// f(x) = j0(x)
static J0_TBL: [ZeroPoly; 15] = [
    // [1, 3.9269908169872414)
    ZeroPoly {
        c: f64::from_bits(0x40033D152E971B40), // 2.404825557695773e0
        k_hi: [
            F64x2::from_bits(0xBC919B7921F03C8E, 0xB90ECE5CFC3B5051), // -6.108765259736730397081979074235e-17
            F64x2::from_bits(0xBFE09CDB36551280, 0xBC784836A8641038), // -5.191474972894667627380839536040e-1
            F64x2::from_bits(0x3FBBA1DEEA029495, 0x3C0A956E29D4C33C), // 1.079387017549201107227730998483e-1
            F64x2::from_bits(0x3FACFAE864368D70, 0xBC339A933C632462), // 5.660177443794622736159855765893e-2
            F64x2::from_bits(0xBF81BB1CBE1A4072, 0x3BF7D893366B8EDA), // -8.657669593304915314466490603787e-3
            F64x2::from_bits(0xBF61F992590D12BD, 0x3C0E9D6F09F949F0), // -2.194200359016155889320141280192e-3
            F64x2::from_bits(0x3F315382BA06CC48, 0x3BCFBF4F1DD63979), // 2.643770367526982220503331224258e-4
            F64x2::from_bits(0x3F06ED3B9F07EB21, 0xBB9F78FC53549F0A), // 4.372919272906659872744439466529e-5
            F64x2::from_bits(0xBED232C77D228B42, 0xBB77B102E90DDABC), // -4.338826288615925756895215664661e-6
            F64x2::from_bits(0xBEA1CCE3028212B1, 0x3B3D1BEA417F7E63), // -5.304914383474427358150217599747e-7
            F64x2::from_bits(0x3E67FF99166CB615, 0xBAF8504EB99AF4F4), // 4.470055863510828286610697242164e-8
        ],
        k_lo: [
            f64::from_bits(0x3E32951BD46D39DC), // 4.3265651676153165e-9
            f64::from_bits(0xBDF5C2C38B87ADF3), // -3.166612627756426e-10
            f64::from_bits(0xBDBBDC4686315F5E), // -2.5338932676961315e-11
            f64::from_bits(0x3D7CD41D38A15CC9), // 1.638714529708127e-12
            f64::from_bits(0x3D3F70ABF0B868A5), // 1.1169775718038169e-13
            f64::from_bits(0xBCFD22F7F4334689), // -6.469623942433395e-15
            f64::from_bits(0xBCBBA868DB6740CF), // -3.838297734026803e-16
            f64::from_bits(0x3C773D4734D07429), // 2.0156938420473082e-17
            f64::from_bits(0x3C33149614C279FB), // 1.0343513727675074e-18
            f64::from_bits(0xBBEF148FD3CA202B), // -5.2652109476417294e-20
        ],
    },
    // [3.9269908169872414, 7.0685834705770345)
    ZeroPoly {
        c: f64::from_bits(0x4016148F5B2C2E45), // 5.520078110286311e0
        k_hi: [
            F64x2::from_bits(0xBC7FBB40985F6E34, 0x3911AC8A875E2310), // -2.752264943262183147205899994392e-17
            F64x2::from_bits(0x3FD5C6E60A097823, 0xBC5EE4B8E4EEB82E), // 3.402648065583681535507795965839e-1
            F64x2::from_bits(0xBF9F8F72E7A848DC, 0xBC3C2F8E2816AB45), // -3.082065142559363577074766193024e-2
            F64x2::from_bits(0xBFAB2150CB41E8C1, 0xBC037ABFE5AA5F5D), // -5.298855286760462842564997671185e-2
            F64x2::from_bits(0x3F72F7FFE90256BA, 0xBC12F7A6C413EC53), // 4.631042145907643493534796694257e-3
            F64x2::from_bits(0x3F627E31FE9A9779, 0x3BCC00D65F7EC1BD), // 2.257440229032900176683762615284e-3
            F64x2::from_bits(0xBF26F641F41956F5, 0xBBC514253A3286B2), // -1.751857289941577294614069824080e-4
            F64x2::from_bits(0xBF0863F481A43012, 0x3B96CFB589F51102), // -4.652109073974298600390249229022e-5
            F64x2::from_bits(0x3ECAD77D748A05E9, 0xBB4649B488C3CDBF), // 3.199786906827403332714514092638e-6
            F64x2::from_bits(0x3EA32E6D99C697E7, 0x3B42ECD76FE9F89B), // 5.716490707921113373523018861679e-7
            F64x2::from_bits(0xBE62DA37E383C00C, 0xBB0876A96D41EEC6), // -3.511535989926808104968006742419e-8
        ],
        k_lo: [
            f64::from_bits(0xBE341D72D927AAFE), // -4.683396242418548e-9
            f64::from_bits(0x3DF1D0433D567C83), // 2.592209201959078e-10
            f64::from_bits(0x3DBE2F3379C73F25), // 2.7452532964884934e-11
            f64::from_bits(0xBD78498FCE8B2670), // -1.3805760413442374e-12
            f64::from_bits(0xBD4104733511F35D), // -1.2091577818273307e-13
            f64::from_bits(0x3CF911949F511742), // 5.566363824716137e-15
            f64::from_bits(0x3CBDE37C47D77BCA), // 4.14787858531783e-16
            f64::from_bits(0xBC744637AEB3DA71), // -1.7585140935848565e-17
            f64::from_bits(0xBC347739075DAE0C), // -1.1094485757618615e-18
            f64::from_bits(0x3BE937DF1E7BFBD9), // 4.272137432526093e-20
        ],
    },
    // [7.0685834705770345, 10.210176124166829)
    ZeroPoly {
        c: f64::from_bits(0x40214EB56CCCDECA), // 8.653727912911013e0
        k_hi: [
            F64x2::from_bits(0xBC96E8EEB22E5818, 0x39363AB447C03C7A), // -7.948465570525161599981923342319e-17
            F64x2::from_bits(0xBFD15F7977A772D4, 0x3C7B4F45A2B20B61), // -2.714522999283819112030025218104e-1
            F64x2::from_bits(0x3F900F7FCF183E18, 0x3C2CCB085B50D3CD), // 1.568412496095392205992130800448e-2
            F64x2::from_bits(0x3FA68B984EC6493B, 0x3C4001845B3D0752), // 4.403377496341184284331792945418e-2
            F64x2::from_bits(0xBF648E63600D841F, 0x3C076B033F995F93), // -2.509302227210699028295558859611e-3
            F64x2::from_bits(0xBF60E0D60385A6F0, 0x3C0D092417C6B5DD), // -2.060335155125929210462430710487e-3
            F64x2::from_bits(0x3F1D796052775AB0, 0x3BBF0790B3AA316B), // 1.124348678936775196079059855854e-4
            F64x2::from_bits(0x3F07800BC5506714, 0x3BAC88858B7D2E29), // 4.482303544524175708147913629565e-5
            F64x2::from_bits(0xBEC3324842F7D46F, 0x3B61F95ECFD42D36), // -2.288391008904150770725978388288e-6
            F64x2::from_bits(0xBEA30E8CC35F07D2, 0x3B3D71ECB9EC66D1), // -5.679379515499026794377271376347e-7
            F64x2::from_bits(0x3E5CEDA479A08F8E, 0xBAE9F71CF565E31C), // 2.694157039835279091624471875359e-8
        ],
        k_lo: [
            f64::from_bits(0x3E34582CB2054D54), // 4.736807196768933e-9
            f64::from_bits(0xBDECADE19EFA38B8), // -2.08668817892119e-10
            f64::from_bits(0xBDBEDAFBDC806B60), // -2.8062827913323888e-11
            f64::from_bits(0x3D743C22EB42109B), // 1.1502213408960854e-12
            f64::from_bits(0x3D4182089E91DA70), // 1.244014244278836e-13
            f64::from_bits(0xBCF56D3B8774B1F6), // -4.757680825385067e-15
            f64::from_bits(0xBCBEDE9A093CA1AE), // -4.2840089675802875e-16
            f64::from_bits(0x3C71AC7F1ACA1986), // 1.5329590429500732e-17
            f64::from_bits(0x3C352EEE14590065), // 1.1483500947358345e-18
            f64::from_bits(0xBBE66C89EA97FDC0), // -3.798769878334432e-20
        ],
    },
    // [10.210176124166829, 13.351768777756622)
    ZeroPoly {
        c: f64::from_bits(0x40279544008272B6), // 1.1791534439014281e1
        k_hi: [
            F64x2::from_bits(0xBC92D8ED368E0843, 0xB93E1AA41646D0F9), // -6.538994895807815285223184197569e-17
            F64x2::from_bits(0x3FCDC13E66AC2E77, 0x3C25DDDAD54782C9), // 2.324598313647247844014354496249e-1
            F64x2::from_bits(0xBF842FF0CDC58451, 0x3C29123484EECAEA), // -9.857064513825759489065457572839e-3
            F64x2::from_bits(0xBFA38D1DD8992E05, 0x3C40B330E6C25F71), // -3.818600911162309776607953779132e-2
            F64x2::from_bits(0x3F5A55E9B346EDB5, 0xBBD5215B10A94200), // 1.607397292089679891537294590848e-3
            F64x2::from_bits(0x3F5E2E16F97D0A81, 0x3BD24DEAB5CBC74D), // 1.842043338868976320196880119053e-3
            F64x2::from_bits(0xBF13DFC3782AF1FD, 0xBBB51517679C01A8), // -7.581358480996527027257227401552e-5
            F64x2::from_bits(0xBF05CE7F496656A9, 0x3B8230F69F5431CA), // -4.159284542308339778274541127569e-5
            F64x2::from_bits(0x3EBBB178DA9C49A0, 0x3B5EAD28F3B2537C), // 1.650646348887002747980389815502e-6
            F64x2::from_bits(0x3EA2346D7493F2BD, 0xBB4B74735FDA44E0), // 5.425451742587530231205635605221e-7
            F64x2::from_bits(0xBE5612F2E7992DA7, 0x3AFC60D5344C05A2), // -2.055803214320044553254829946550e-8
        ],
        k_lo: [
            f64::from_bits(0xBE33D79882D33394), // -4.61986541834555e-9
            f64::from_bits(0x3DE6DBC1126E2688), // 1.6631685717452698e-10
            f64::from_bits(0x3DBE88C4C7AEAC90), // 2.7770740984715337e-11
            f64::from_bits(0xBD70AD87F3818E16), // -9.480263372902583e-13
            f64::from_bits(0xBD4180D18CE8249F), // -1.2436769825977725e-13
            f64::from_bits(0x3CF22087C16E3EF6), // 4.025018422091023e-15
            f64::from_bits(0x3CBF156ACAFA7CF9), // 4.3137244859631486e-16
            f64::from_bits(0xBC6E876D657BB8A4), // -1.3239848890360323e-17
            f64::from_bits(0xBC356FC7187A3EF5), // -1.1620821334291275e-18
            f64::from_bits(0x3BE3B9F0B3BECFAC), // 3.341770130647886e-20
        ],
    },
    // [13.351768777756622, 16.493361431346415)
    ZeroPoly {
        c: f64::from_bits(0x402DDCA13EF271D2), // 1.4930917708487787e1
        k_hi: [
            F64x2::from_bits(0xBCA50BE2EF09843E, 0xB93ED6A5C02E698D), // -1.460389551237032993995631139785e-16
            F64x2::from_bits(0xBFCA701D0F967500, 0x3C68CDD47D8A73B0), // -2.065464330779960170480261480711e-1
            F64x2::from_bits(0x3F7C54B930FEF8E5, 0x3C00B335D01A886B), // 6.916736034268765607378964025673e-3
            F64x2::from_bits(0x3FA17798AA09F11F, 0xBC332079AF955C5D), // 3.411557269734743249153413195386e-2
            F64x2::from_bits(0xBF52A2151407DD25, 0x3BFAF577967C8D0A), // -1.137276252959848071067528494756e-3
            F64x2::from_bits(0xBF5B541F829BFB40, 0x3BDAFBBDEB62568A), // -1.668005725548218886504031580837e-3
            F64x2::from_bits(0x3F0CC0BDA19FFE4E, 0x3B98B5AE2FA073DA), // 5.484179219194790483590089981645e-5
            F64x2::from_bits(0x3F041F3B0662F476, 0xBB99B0A5E1ED6ECC), // 3.837965850681742562261212152807e-5
            F64x2::from_bits(0xBEB4B230B114B070, 0x3B498F2DE2739837), // -1.233581034941327952101776854567e-6
            F64x2::from_bits(0xBEA1223E2BD04F62, 0x3B34924CA1CD3CF8), // -5.106258736974767416181869098856e-7
            F64x2::from_bits(0x3E5119639689D1E1, 0xBAF9820F16A9151D), // 1.592484846981482983235394911568e-8
        ],
        k_lo: [
            f64::from_bits(0x3E32FFB8DEC7E98A), // 4.423529525744808e-9
            f64::from_bits(0xBDE24C47F1828D22), // -1.3313527429373914e-10
            f64::from_bits(0xBDBDA527431704D0), // -2.6962088979093204e-11
            f64::from_bits(0x3D6B75EE83610BA6), // 7.804792026206686e-13
            f64::from_bits(0x3D412DCADA5D16D4), // 1.220632593463989e-13
            f64::from_bits(0xBCEE8DB4B4FE7F25), // -3.39212420598333e-15
            f64::from_bits(0xBCBEC3AAF49FE23E), // -4.2694080666510743e-16
            f64::from_bits(0x3C6A3BC423CF023F), // 1.1376950426457076e-17
            f64::from_bits(0x3C3559C3E19A2DA7), // 1.15742079336585e-18
            f64::from_bits(0xBBE13D82EE29D30D), // -2.9206168815166026e-20
        ],
    },
    // [16.493361431346415, 19.634954084936208)
    ZeroPoly {
        c: f64::from_bits(0x403212313F8A19F6), // 1.8071063967910924e1
        k_hi: [
            F64x2::from_bits(0x3CAA2122AF76659E, 0x394C77A24A5F052B), // 1.813093807534984955199404354603e-16
            F64x2::from_bits(0x3FC8077F56C9B782, 0xBC608EA6507D275D), // 1.877288030404394220477625012460e-1
            F64x2::from_bits(0xBF75467EB535DF12, 0x3BD75F926F1F133A), // -5.194182350684932774806463527426e-3
            F64x2::from_bits(0xBF9FD7C3AD6F59DF, 0x3C22CED7B84459C1), // -3.109651323378581235326125997359e-2
            F64x2::from_bits(0x3F4C1B47C809C62F, 0x3BECDBF96D712574), // 8.577442641366002472386271063808e-4
            F64x2::from_bits(0x3F59166C7D3EAA7B, 0xBBAA31E4BCD94AD8), // 1.531225153466180289795595112748e-3
            F64x2::from_bits(0xBF05F0152478F618, 0xBB92C8FB6447CC2E), // -4.184307595920757891785079973853e-5
            F64x2::from_bits(0xBF02AA939FD9FC73, 0x3BA6C9F0905FE97B), // -3.560317052819167408349018940281e-5
            F64x2::from_bits(0x3EB0129708EC235F, 0x3B51159CF1F37F76), // 9.580026333907722458944278182168e-7
            F64x2::from_bits(0x3EA01716DC9F1517, 0xBB440FFDEE5645E3), // 4.795251069492244098890988776489e-7
            F64x2::from_bits(0xBE4B227EEF3FDACE, 0xBAEBD5F9A018CBD8), // -1.263560232138229073852162323517e-8
        ],
        k_lo: [
            f64::from_bits(0xBE32101C49C384C6), // -4.205604000658624e-9
            f64::from_bits(0x3DDDB1E8B73E50FB), // 1.0802962501191958e-10
            f64::from_bits(0x3DBC80AAA10C1DAF), // 2.5922966954848924e-11
            f64::from_bits(0xBD66C2B152E7F0E8), // -6.468928361086886e-13
            f64::from_bits(0xBD40AD5AA7FC9D71), // -1.1849838126783977e-13
            f64::from_bits(0x3CE9D01E763BF8EE), // 2.8658147866992686e-15
            f64::from_bits(0x3CBE1BA2BC64C843), // 4.178317678372532e-16
            f64::from_bits(0xBC6689E44852D0F2), // -9.77457679237609e-18
            f64::from_bits(0xBC3508003571BA63), // -1.1401065196919745e-18
            f64::from_bits(0x3BDE11AE7360A139), // 2.546949142164612e-20
        ],
    },
    // [19.634954084936208, 22.776546738526)
    ZeroPoly {
        c: f64::from_bits(0x4035362DD173F792), // 2.121163662987926e1
        k_hi: [
            F64x2::from_bits(0x3C98B4B912DA2218, 0x391C38E61A6B7DE1), // 8.571597945195108970368302472838e-17
            F64x2::from_bits(0xBFC62D93AA9D05BB, 0xBC6D962D0BF4C4ED), // -1.732658942292298737767046015338e-1
            F64x2::from_bits(0x3F70BA9CE88929C0, 0x3C18B5A1E4C95422), // 4.084217951979274620765878028220e-3
            F64x2::from_bits(0x3F9D7073DAEBB02C, 0x3C21D2BFE5CC2286), // 2.874928497014665788615607254604e-2
            F64x2::from_bits(0xBF462813C7F58712, 0xBBE83C5E53B74880), // -6.761643016199601006126694187337e-4
            F64x2::from_bits(0xBF574A948D05638E, 0xBBF02770AFD46120), // -1.421589917376151239373536915173e-3
            F64x2::from_bits(0x3F01695765F13BAD, 0xBB8636EDDCE9960E), // 3.320978134323886763309038876506e-5
            F64x2::from_bits(0x3F0170AB5EEB0E28, 0x3BA5B6EC883E9A9B), // 3.326437933062433983831323874314e-5
            F64x2::from_bits(0xBEA9C8E418A201B2, 0xBB4D75909568A3C2), // -7.684448502104828475223381203216e-7
            F64x2::from_bits(0xBE9E4D893DCE022A, 0x3B2AAA0428093BBD), // -4.515480307929213919464522718427e-7
            F64x2::from_bits(0x3E460F73A6CC0D04, 0xBAE02609BF3E5836), // 1.027265491551836406150101930213e-8
        ],
        k_lo: [
            f64::from_bits(0x3E312782F0E0059A), // 3.994056431076496e-9
            f64::from_bits(0xBDD881D583CF3FDF), // -8.91565440853185e-11
            f64::from_bits(0xBDBB4C5209AB273B), // -2.4827501690192228e-11
            f64::from_bits(0x3D63132D8488AF9D), // 5.421416430847788e-13
            f64::from_bits(0x3D401931A8C872CE), // 1.1438611118503365e-13
            f64::from_bits(0xBCE5F49C48F10297), // -2.4375512407129556e-15
            f64::from_bits(0xBCBD44A0E3FC93F0), // -4.0617620371708233e-16
            f64::from_bits(0x3C636F7DBE1C4CA4), // 8.428809920035235e-18
            f64::from_bits(0x3C34913DE76298E1), // 1.1149583224865658e-18
            f64::from_bits(0xBBDA41CDC61BDC06), // -2.2240583152714882e-20
        ],
    },
    // [22.776546738526, 25.918139392115794)
    ZeroPoly {
        c: f64::from_bits(0x40385A3B930156DD), // 2.4352471530749302e1
        k_hi: [
            F64x2::from_bits(0xBCA55E059345B430, 0xB94E638E16E676F8), // -1.482652373933724548746063372456e-16
            F64x2::from_bits(0x3FC4B2A2EBF61ECE, 0xBC6053B468D823F4), // 1.617015506892500105421159052085e-1
            F64x2::from_bits(0xBF6B3297FDAE7857, 0xBC0F1ADCFB268863), // -3.320023400603704571115108052184e-3
            F64x2::from_bits(0xBF9B8105D59B114C, 0xBC35D33B4134DF38), // -2.685937038656165239478549250790e-2
            F64x2::from_bits(0x3F420A3F8C12A1C7, 0xBBE14E48C9CE3A8F), // 5.505380905975125829101531113674e-4
            F64x2::from_bits(0x3F55D18D69DE6CFB, 0x3BE166684BC183AA), // 1.331699465913881497134238065098e-3
            F64x2::from_bits(0xBEFC79DB4B341E50, 0xBB8A7B25EBBBBB63), // -2.715683212975293056244115869127e-5
            F64x2::from_bits(0xBF00679C92C3038D, 0xBB9C9F856414C100), // -3.128954480891626917022949663816e-5
            F64x2::from_bits(0x3EA53AC4ECD4BCE9, 0xBB337179F6900DDA), // 6.326904099847175387552159342300e-7
            F64x2::from_bits(0x3E9CA75078409025, 0x3B1B20031608446F), // 4.269714894997273877951461273604e-7
            F64x2::from_bits(0xBE4253415F08E804, 0x3AB264D8A04AA5B9), // -8.533343782557164395867304270287e-9
        ],
        k_lo: [
            f64::from_bits(0xBE305141504AC3DC), // -3.799191409970691e-9
            f64::from_bits(0x3DD49199409D8462), // 7.482865729288693e-11
            f64::from_bits(0x3DBA20DE854C4F68), // 2.376363717881484e-11
            f64::from_bits(0xBD602FFF1FFD4C09), // -4.600760419160176e-13
            f64::from_bits(0xBD3F01C86EF45F9D), // -1.1015886734789526e-13
            f64::from_bits(0x3CE2D780ECD9C653), // 2.091861238939873e-15
            f64::from_bits(0x3CBC59AF1AEF9587), // 3.9343983825864154e-16
            f64::from_bits(0xBC60DC94DA49B90B), // -7.312573425156409e-18
            f64::from_bits(0xBC3406A51E4C4711), // -1.0856093046054439e-18
            f64::from_bits(0x3BD703EB8D72C097), // 1.9494728399322973e-20
        ],
    },
    // [25.918139392115794, 29.059732045705587)
    ZeroPoly {
        c: f64::from_bits(0x403B7E54A5FD5F11), // 2.7493479132040253e1
        k_hi: [
            F64x2::from_bits(0x3CB1C17ABE35EAAE, 0x3947A7A9046DA289), // 2.464109355407528432594490269251e-16
            F64x2::from_bits(0xBFC37AAC8C1AEABB, 0xBC500CB2021E9216), // -1.521812137705945457252074464730e-1
            F64x2::from_bits(0x3F66AC0D2E2F2E6C, 0x3BF095847F089DAE), // 2.767587416632902017744070344947e-3
            F64x2::from_bits(0x3F99E74E754EA71F, 0x3C38A594E4C3660F), // 2.529642668711240336249755558451e-2
            F64x2::from_bits(0xBF3E1C0589E32AF5, 0xBBC27EB5ABAE3B7D), // -4.594338914897655112244880497516e-4
            F64x2::from_bits(0xBF5496158DC5F7FF, 0xBBF427812E6FF406), // -1.256485985407618095825541489622e-3
            F64x2::from_bits(0x3EF7D55405348C44, 0xBB896210116B42F0), // 2.272921876768039463540204747060e-5
            F64x2::from_bits(0x3EFF0B30F45061E1, 0x3B91416E3CC32A11), // 2.960559438145046042647299605176e-5
            F64x2::from_bits(0xBEA1D9E3629B97B7, 0x3B45FD581C50864B), // -5.320050082385274136225145091873e-7
            F64x2::from_bits(0xBE9B35EF28E9AF2A, 0xBB3E52CB3DB6047A), // -4.054707369086232924021631484808e-7
            F64x2::from_bits(0x3E3F0200E2640F1C, 0x3AD163D1D193FAD8), // 7.219572084477759540433136396198e-9
        ],
        k_lo: [
            f64::from_bits(0x3E2F1FB225C57254), // 3.623288598120748e-9
            f64::from_bits(0xBDD187928A1C59F2), // -6.377223968318564e-11
            f64::from_bits(0xBDB908BA748160C9), // -2.2768376838055335e-11
            f64::from_bits(0x3D5BD0A9402B1B20), // 3.9527609715272e-13
            f64::from_bits(0x3D3DD9813AC4AD11), // 1.0604718218818714e-13
            f64::from_bits(0xBCE0535D79EB8112), // -1.812510706453432e-15
            f64::from_bits(0xBCBB6B79C14D6820), // -3.805265350764151e-16
            f64::from_bits(0x3C5D784C990BCF47), // 6.390269995067546e-18
            f64::from_bits(0x3C3373EF48F7F26C), // 1.0545421930352684e-18
            f64::from_bits(0xBBD447A75000825B), // -1.7177740702530204e-20
        ],
    },
    // [29.059732045705587, 32.20132469929538)
    ZeroPoly {
        c: f64::from_bits(0x403EA27591CBBED2), // 3.0634606468431976e1
        k_hi: [
            F64x2::from_bits(0x3C966608AC164DBD, 0x393546DAC453BC2E), // 7.771064981615525857208941779875e-17
            F64x2::from_bits(0x3FC27407DFADEE6D, 0xBC51D9F7B52D2549), // 1.441659776863732051059756724005e-1
            F64x2::from_bits(0xBF6346950BFD924A, 0xBBF037445FB8D523), // -2.352992159943921807992950335790e-3
            F64x2::from_bits(0xBF988D48D1D4EB7A, 0xBC328D0F0F168B20), // -2.397645741187257834671329782060e-2
            F64x2::from_bits(0x3F399E6923ADA95E, 0xBBDF3143533C88A8), // 3.909117401038798190147225471026e-4
            F64x2::from_bits(0x3F538984B76CDF4A, 0xBBEB9F2537154F40), // 1.192454915964052298813157455713e-3
            F64x2::from_bits(0xBEF452194B75F3E4, 0xBB6B03F747B32EC2), // -1.937932821992676103640255167306e-5
            F64x2::from_bits(0xBEFD855D7B8F01FC, 0x3B9AF7D5D88589C5), // -2.815337913708754081446495897944e-5
            F64x2::from_bits(0x3E9E8A88601FF0FA, 0xBB221872090988C4), // 4.550985012172954499030771743399e-7
            F64x2::from_bits(0x3E99F1BD69B140A8, 0xBB1190384D643F8B), // 3.866001436246982532116991930449e-7
            F64x2::from_bits(0xBE3AA483FCDF7954, 0x3AA2452EC5439558), // -6.203222780900367440967513192717e-9
        ],
        k_lo: [
            f64::from_bits(0xBE2DC4C0236FD5B5), // -3.4655161201605962e-9
            f64::from_bits(0x3DCE48DC1ACDD465), // 5.508738200885631e-11
            f64::from_bits(0x3DB807111F0897B2), // 2.1852979443027527e-11
            f64::from_bits(0xBD582C0E31692096), // -3.4350608143755783e-13
            f64::from_bits(0xBD3CC15084660DE7), // -1.0215876089851676e-13
            f64::from_bits(0x3CDC8F158BD0BD3D), // 1.5853386669105086e-15
            f64::from_bits(0x3CBA835E1FD0055C), // 3.679439388219673e-16
            f64::from_bits(0xBC59F2F68D9EAE80), // -5.6268086117794776e-18
            f64::from_bits(0xBC32E080B4559726), // -1.0233222620932728e-18
            f64::from_bits(0x3BD1F9E344EA2395), // 1.5226369378959504e-20
        ],
    },
    // [32.20132469929538, 35.34291735288517)
    ZeroPoly {
        c: f64::from_bits(0x4040E34E13A66FE6), // 3.377582021357357e1
        k_hi: [
            F64x2::from_bits(0x3CACC62C8409DAA7, 0xB94DD62BD7ACEE57), // 1.996605433542561316731093339327e-16
            F64x2::from_bits(0xBFC192F23CE3E051, 0x3C61BDE0BCC76480), // -1.372969434085029783672519646809e-1
            F64x2::from_bits(0x3F60A668185C00CB, 0x3C0157107B0381BE), // 2.032473860595106932248101884648e-3
            F64x2::from_bits(0x3F9764141D65208A, 0xBBFF1E42A3663D8B), // 2.284270696080218505112344391203e-2
            F64x2::from_bits(0xBF3624437A2FE6D1, 0xBBC9FD1AC40FD3BD), // -3.378548373685732630811432876283e-4
            F64x2::from_bits(0xBF52A184BE0D9892, 0x3BF936743C8B44B7), // -1.137141829724384532653749069309e-3
            F64x2::from_bits(0x3EF196DE0EEEF13E, 0xBB8CB6B05E339D46), // 1.677448829771995805692843579900e-5
            F64x2::from_bits(0x3EFC317F85411269, 0x3B82DD388180EF28), // 2.688727575185721295128631391656e-5
            F64x2::from_bits(0xBE9A8019EF771FAD, 0x3B133D2E27FADEF8), // -3.948866686667847499459489461449e-7
            F64x2::from_bits(0xBE98D38497BEBDAB, 0xBB3E9E76C1D68C39), // -3.699398332069506807713162846455e-7
            F64x2::from_bits(0x3E37318410807B69, 0x3ADCF4D453AEF279), // 5.400139231641108934907430104980e-9
        ],
        k_lo: [
            f64::from_bits(0x3E2C8D9D45B65C69), // 3.3240277593795356e-9
            f64::from_bits(0xBDCA782C1A5AA4E0), // -4.814759987774637e-11
            f64::from_bits(0xBDB71BBB6603E739), // -2.101690208668264e-11
            f64::from_bits(0x3D5539500292D5DA), // 3.016094340540869e-13
            f64::from_bits(0x3D3BBC1A6604E87B), // 9.853372450170694e-14
            f64::from_bits(0xBCD933941F7174B8), // -1.398963107853771e-15
            f64::from_bits(0xBCB9A61D1735756B), // -3.5594973341374705e-16
            f64::from_bits(0x3C57050169DC8BEF), // 4.991569835119993e-18
            f64::from_bits(0x3C3250AFA8B4F69B), // 9.928679159105918e-19
            f64::from_bits(0xBBD0083BD4403B56), // -1.3579770208729947e-20
        ],
    },
    // [35.34291735288517, 38.48451000647496)
    ZeroPoly {
        c: f64::from_bits(0x404275637A9619EC), // 3.6917098353664045e1
        k_hi: [
            F64x2::from_bits(0x3CA18E7B1F5C77AE, 0x3945E73CE80C9B6A), // 1.218231501809992623682780729341e-16
            F64x2::from_bits(0x3FC0CF3ED059C573, 0x3C32E2CC445527D1), // 1.313246266686679274672983284365e-1
            F64x2::from_bits(0xBF5D242AA5299435, 0x3BB3A868BEFE12D3), // -1.778642316503165452316118635932e-3
            F64x2::from_bits(0xBF96613D93B0180B, 0x3C352A2E78388318), // -2.185531820216372416339366344996e-2
            F64x2::from_bits(0x3F33627F261F5174, 0xBBD8D39B95028D3D), // 2.957878515392954822630510280296e-4
            F64x2::from_bits(0x3F51D69CA0D88394, 0xBBEB0C9BFC67ADC1), // 1.088765103652832150178522281839e-3
            F64x2::from_bits(0xBEEED574AFAB715D, 0x3B827947C0939D03), // -1.470270716679631616418092709628e-5
            F64x2::from_bits(0xBEFB06384E48EE17, 0xBB8B317818AE5FFE), // -2.577237764199270346617501038449e-5
            F64x2::from_bits(0x3E97452C652355FB, 0x3B352B451BB81984), // 3.467531303318611502385063672925e-7
            F64x2::from_bits(0x3E97D51F133B3CEC, 0x3AB75BFD8D52B2B9), // 3.551320048979163815915585256060e-7
            F64x2::from_bits(0xBE346AC67A7D5DD9, 0x3AD2FA8D836EDFE1), // -4.753724448659974370848489538759e-9
        ],
        k_lo: [
            f64::from_bits(0xBE2B75D66A6D5C7C), // -3.196800008231585e-9
            f64::from_bits(0x3DC75F5A5D4EA4D0), // 4.251427999871935e-11
            f64::from_bits(0x3DB64531BBC1B1A4), // 2.0254710872750832e-11
            f64::from_bits(0xBD52CEABF011A2BA), // -2.67267965087833e-13
            f64::from_bits(0xBD3ACA415271D0E7), // -9.517740990137395e-14
            f64::from_bits(0x3CD66B86AB477E8B), // 1.2445613227093753e-15
            f64::from_bits(0x3CB8D5C0C7A30368), // 3.446544832416521e-16
            f64::from_bits(0xBC5490B3BDADD5B4), // -4.459376147209385e-18
            f64::from_bits(0xBC31C6D45CDC9D50), // -9.63675639844207e-19
            f64::from_bits(0x3BCCC50F63460C62), // 1.2184469864152782e-20
        ],
    },
    // [38.48451000647496, 41.62610266006476)
    ZeroPoly {
        c: f64::from_bits(0x4044077A7ED6293A), // 4.005842576462824e1
        k_hi: [
            F64x2::from_bits(0xBC8376CC57901722, 0x38F546E0D1599D40), // -3.376475403459293569958647874660e-17
            F64x2::from_bits(0xBFC0230B9797A7B3, 0x3C6845CE1252F59C), // -1.260694971272734060879617525265e-1
            F64x2::from_bits(0x3F59C8083B2B7588, 0xBBE2E9CBB2AE9A33), // 1.573570287909250057945707847041e-3
            F64x2::from_bits(0x3F957D3203BEFD8F, 0x3C3E7E3CF8AC7C0A), // 2.098539493431200538156230494306e-2
            F64x2::from_bits(0xBF3127CBA22892F8, 0xBBC1BCBB2EAF69C6), // -2.617714073116918129149518787598e-4
            F64x2::from_bits(0xBF51234471455A6B, 0xBBFC0D7B613027D6), // -1.046006048530802017334990926160e-3
            F64x2::from_bits(0x3EEB4FE26EC3E4A1, 0x3B7D1296191AF1F7), // 1.302339975283096475186665912796e-5
            F64x2::from_bits(0x3EF9FC5254F1082A, 0xBB8B9884C5101903), // 2.478182915078439298538887861431e-5
            F64x2::from_bits(0xBE94A44C6506C30E, 0xBADF4B2CE38E7913), // -3.075866503955760211997700459988e-7
            F64x2::from_bits(0xBE96F17DD184832A, 0x3B261435CE8D37EC), // -3.418822003432587970203805455999e-7
            F64x2::from_bits(0x3E3225640A69F5B8, 0xBAD30A9D573A4DC6), // 4.224958305625103922048648172434e-9
        ],
        k_lo: [
            f64::from_bits(0x3E2A7943503DCDF7), // 3.081942369652535e-9
            f64::from_bits(0xBDC4D296EC8EB81A), // -3.7876116799970865e-11
            f64::from_bits(0xBDB58176F7712974), // -1.955933979461788e-11
            f64::from_bits(0x3D50CD719D5882CA), // 2.387780976917688e-13
            f64::from_bits(0x3D39EAF2CED76A23), // 9.207840696371174e-14
            f64::from_bits(0xBCD417D072CD3B61), // -1.1153869173185388e-15
            f64::from_bits(0xBCB812CDF0259729), // -3.3408629842613907e-16
            f64::from_bits(0x3C527E4DF6A0E8A5), // 4.01011193152982e-18
            f64::from_bits(0x3C31440AD326E347), // 9.35980360790731e-19
            f64::from_bits(0xBBC9F67163AC5D27), // -1.0995617463547102e-20
        ],
    },
    // [41.62610266006476, 44.76769531365455)
    ZeroPoly {
        c: f64::from_bits(0x40459992C65D0D8D), // 4.319979171317673e1
        k_hi: [
            F64x2::from_bits(0xBC9DF12B66B26573, 0x39017C7E109B326D), // -1.038824182895261314033313223307e-16
            F64x2::from_bits(0x3FBF13FAF32C8E0B, 0xBC5CCD3634EB6D51), // 1.213986247717501541772809486895e-1
            F64x2::from_bits(0xBF570558DDDB7C57, 0x3BC50F3E38B4C775), // -1.405083450144469763458589887717e-3
            F64x2::from_bits(0xBF94B24D7A933973, 0x3C3EDF07430C8D5E), // -2.021142063750276148550107373880e-2
            F64x2::from_bits(0x3F2EA52A21487972, 0xBBC1DD42AD8391ED), // 2.338041237103894180875805711052e-4
            F64x2::from_bits(0x3F50834D8F3FDD5B, 0xBBEC6FD7CCF6CF72), // 1.007867566839090989246068916538e-3
            F64x2::from_bits(0xBEE86941A4B43B93, 0x3B8FD4228DEBD142), // -1.164014715595936879693474324345e-5
            F64x2::from_bits(0xBEF90E32CEF3E8C1, 0xBB44FD6C7A4F0716), // -2.389475133399914133007253004026e-5
            F64x2::from_bits(0x3E92785F6385B120, 0xBAE3D42AB792E529), // 2.752275096915414941950234908572e-7
            F64x2::from_bits(0x3E962482BF9D02A5, 0xBB3E045A26EE35DB), // 3.299507508221580337606250379411e-7
            F64x2::from_bits(0xBE30431253861DAB, 0xBAAED8CEDA5C4192), // -3.786291551453743754298031641348e-9
        ],
        k_lo: [
            f64::from_bits(0xBE29943C301FBB3D), // -2.9777925689419198e-9
            f64::from_bits(0x3DC2B24579C5F55B), // 3.400850367088081e-11
            f64::from_bits(0x3DB4CE82A200D789), // 1.8923565947922367e-11
            f64::from_bits(0xBD4E3DCDA85EC203), // -2.1487820829265067e-13
            f64::from_bits(0xBD391CD773694DD4), // -8.921809964128822e-14
            f64::from_bits(0x3CD220B3247E6302), // 1.0062913557116253e-15
            f64::from_bits(0x3CB75CF943404693), // 3.2422923300291937e-16
            f64::from_bits(0xBC50BB12B38597EF), // -3.6279039903597955e-18
            f64::from_bits(0xBC30C8B24859D2B5), // -9.098608572784269e-19
            f64::from_bits(0x3BC78E6B833BE6BC), // 9.97649297149446e-21
        ],
    },
    // [44.76769531365455, 47.909287967244346)
    ZeroPoly {
        c: f64::from_bits(0x40472BAC0F810810), // 4.6341188371661815e1
        k_hi: [
            F64x2::from_bits(0xBC9EAB3B2BF5813D, 0xB924FAA01504D293), // -1.064040232267980066021507031278e-16
            F64x2::from_bits(0xBFBE018D99F5DA1B, 0xBC4F5DADCE6D1AD4), // -1.172111988906653807374962663807e-1
            F64x2::from_bits(0x3F54B85897B3635A, 0x3BEFAC6A11AF332A), // 1.264654651825311700247173704141e-3
            F64x2::from_bits(0x3F93FC442153435C, 0x3C3EFC23E2FCD299), // 1.951700642968333566977320669954e-2
            F64x2::from_bits(0xBF2B9694D7148787, 0x3BCFCB61A1EC930B), // -2.104813279766826569992123155831e-4
            F64x2::from_bits(0xBF4FE6FDC644DDDE, 0x3BEAD30CE5A637B9), // -9.735812314265975139856325587371e-4
            F64x2::from_bits(0x3EE5FD096E452450, 0xBB82D591303B19B2), // 1.048489816270390853211077202833e-5
            F64x2::from_bits(0x3EF83770C9A84459, 0x3B88156654CBAA01), // 2.309471583756094979005591709244e-5
            F64x2::from_bits(0xBE90A6F6F7E05E43, 0x3B1EA956E16E4270), // -2.481372048958603860172743572402e-7
            F64x2::from_bits(0xBE956AD4A35E88A5, 0xBB1513D5EF190018), // -3.191427454478845992958929948442e-7
            F64x2::from_bits(0x3E2D5BEA474BA189, 0xBAC1EBC74F12179E), // 3.417842504123646847521669704617e-9
        ],
        k_lo: [
            f64::from_bits(0x3E28C39F68B46BA3), // 2.882926624542033e-9
            f64::from_bits(0xBDC0E7340882B451), // -3.0746628388988736e-11
            f64::from_bits(0xBDB42A6DE03910CA), // -1.8340632843987676e-11
            f64::from_bits(0x3D4B640F5EBAD435), // 1.9462376263967973e-13
            f64::from_bits(0x3D385E6D70257148), // 8.657557299488086e-14
            f64::from_bits(0xBCD073FD27093D6D), // -9.133294981625462e-16
            f64::from_bits(0xBCB6B38FC00708A4), // -3.1504538148493254e-16
            f64::from_bits(0x3C4E7089ECAA621D), // 3.300268539405926e-18
            f64::from_bits(0x3C3054BD51AD096B), // 8.853060309256634e-19
            f64::from_bits(0xBBC57B0717518968), // -9.097378201939717e-21
        ],
    },
];

// GENERATE: bessel_table j1 15 0 2.375 11 10
// J1_TBL[i] has the coefficients `K0, ..., K20` such that
// f(c + d) ~= K0 + K1 * d + ... + K20 * d^20, where `c` is very close to the
// (i + 1)-th positive zero of j1, for the values of the interval of that
// zero, [(i + 0.75) * π, (i + 1.75) * π) (the first one starts at 2.375)
// f(x) = j1(x)
static J1_TBL: [ZeroPoly; 15] = [
    // [2.375, 5.497787143782138)
    ZeroPoly {
        c: f64::from_bits(0x400EA75575AF6F09), // 3.8317059702075125e0
        k_hi: [
            F64x2::from_bits(0xBC91B9C1C3FB286F, 0xB92A19826DC75C7B), // -6.149807356994906091388455188718e-17
            F64x2::from_bits(0xBFD9C6CF582CBF7F, 0x3C7BF01893C9FCB2), // -4.027593957025529560462108644080e-1
            F64x2::from_bits(0x3FAAE8A39F51AD07, 0x3C35E24A68098BE7), // 5.255614585697725973218135252102e-2
            F64x2::from_bits(0x3FAB589D1DA13904, 0x3C36B632F6B9A11D), // 5.341044413272480596451881727478e-2
            F64x2::from_bits(0xBF7537544C331DA9, 0x3C10603591055561), // -5.179719245638571924306032540536e-3
            F64x2::from_bits(0xBF624B3409959064, 0x3C0962CDE05D596F), // -2.233125339147478025166241936852e-3
            F64x2::from_bits(0x3F26E4C2D535422A, 0x3B7EE0A60DF87B56), // 1.746642907201207137276742627373e-4
            F64x2::from_bits(0x3F083A06E30C4103, 0x3BA58EB5B66F75D9), // 4.620870129745833539572529111934e-5
            F64x2::from_bits(0xBEC9799D4C9F28DC, 0x3B60E3ACE7FE73FD), // -3.036863380300169970772906379021e-6
            F64x2::from_bits(0xBEA33825CD2E272D, 0xBB022FCA1C3D4D8A), // -5.727805735139051464718142657206e-7
            F64x2::from_bits(0x3E61706923404C12, 0xBABFDBECA0B54207), // 3.248286298287003547592796243992e-8
        ],
        k_lo: [
            f64::from_bits(0x3E34569B22AAB4C8), // 4.735380569148864e-9
            f64::from_bits(0xBDF03B9E974DCDD8), // -2.3621961964946836e-10
            f64::from_bits(0xBDBEC6230A80A480), // -2.798876465350583e-11
            f64::from_bits(0x3D75EC8593266901), // 1.2462301805141555e-12
            f64::from_bits(0x3D417A3E7F064878), // 1.241852211547798e-13
            f64::from_bits(0xBCF67CF95B26E68B), // -4.9933790150628e-15
            f64::from_bits(0xBCBEE5FB44B9C62C), // -4.288009572753831e-16
            f64::from_bits(0x3C722F8B54F71DC4), // 1.577359752016034e-17
            f64::from_bits(0x3C3563F6EBC2371D), // 1.1595805947623061e-18
            f64::from_bits(0xBBEA53F41794E61B), // -4.460127081726629e-20
        ],
    },
    // [5.497787143782138, 8.63937979737193)
    ZeroPoly {
        c: f64::from_bits(0x401C0FF5F3B47250), // 7.015586669815619e0
        k_hi: [
            F64x2::from_bits(0x3C8049770CE74C2E, 0x390F3497C929687A), // 2.825339409478929375555869776139e-17
            F64x2::from_bits(0x3FD33518B3874E8A, 0x3C764D8939952206), // 3.001157525261325591793517555692e-1
            F64x2::from_bits(0xBF95E70DC60362C3, 0x3BF6E6457CE5EAA9), // -2.138921280934159486210730381341e-2
            F64x2::from_bits(0xBFA80C83BDEEE5B0, 0xBC390EEA469CA1D1), // -4.697047894974149120755675632590e-2
            F64x2::from_bits(0x3F69A4B292E3DE45, 0xBC0849A84747EC13), // 3.130291726048092032380073780837e-3
            F64x2::from_bits(0x3F613FBC7D698217, 0xBBFA003141E29CA9), // 2.105587143248276343928674536880e-3
            F64x2::from_bits(0xBF207358BBDBFF95, 0xBBB90F42CA639117), // -1.255079095514559567037852650977e-4
            F64x2::from_bits(0xBF0796A751F89034, 0x3BA02917AF6F53FA), // -4.499147530293047907615271095935e-5
            F64x2::from_bits(0x3EC4255B015AE161, 0xBB4A1D01CAA91C28), // 2.401580796309598775263987607391e-6
            F64x2::from_bits(0x3EA3026E0CE9677D, 0xBB457366FF55AD57), // 5.665270011543116905919765588791e-7
            F64x2::from_bits(0xBE5D48DCDAEB87C7, 0xBAFFEB432524C6EA), // -2.727342767850030478130756298649e-8
        ],
        k_lo: [
            f64::from_bits(0xBE344639D7DFEB6B), // -4.720483002744274e-9
            f64::from_bits(0x3DEC62CCB5E6EA0D), // 2.0653486815843913e-10
            f64::from_bits(0x3DBECAE9209C34DD), // 2.8005724508666257e-11
            f64::from_bits(0xBD73BB6901A277A2), // -1.1216383779971907e-12
            f64::from_bits(0xBD4183E9A6E64838), // -1.2445357808018774e-13
            f64::from_bits(0x3CF4AE66484460FB), // 4.5921595868002564e-15
            f64::from_bits(0x3CBEF8C87853DE61), // 4.2982019235328446e-16
            f64::from_bits(0xBC70FA875F84CD46), // -1.472661240339588e-17
            f64::from_bits(0xBC3557ED7BC2D4C9), // -1.157031689411286e-18
            f64::from_bits(0x3BE778BBA972CA52), // 3.9762459908369226e-20
        ],
    },
    // [8.63937979737193, 11.780972450961723)
    ZeroPoly {
        c: f64::from_bits(0x402458D0D0BDFC29), // 1.0173468135062722e1
        k_hi: [
            F64x2::from_bits(0x3CA0212F4E592523, 0xB93E073372D2A4BD), // 1.119217779774468185607881890370e-16
            F64x2::from_bits(0xBFCFF654544EBCD2, 0x3C6E78A97B88077B), // -2.497048770578432087509472132198e-1
            F64x2::from_bits(0x3F89223FF2C0783C, 0xBC2929C8D64CB33A), // 1.227235755510146867981852730732e-2
            F64x2::from_bits(0x3FA4B0C5D5DA6789, 0x3C46B77253BB3BA4), // 4.041116939079276545278764995941e-2
            F64x2::from_bits(0xBF5F91A9EE0D2883, 0xBBF0A8F9C4F838B3), // -1.926818797260782686087101488538e-3
            F64x2::from_bits(0xBF5F51C2489B9E6F, 0xBBFB0BD866AFCCEF), // -1.911582689382653479610772027561e-3
            F64x2::from_bits(0x3F16B4C9CA0F7708, 0xBBB064EFC3439514), // 8.661729454170427003810559539537e-5
            F64x2::from_bits(0x3F063C5475439C90, 0x3B9821097E3FBA32), // 4.241116281083437972922213872503e-5
            F64x2::from_bits(0xBEBE3725DAF69C06, 0xBB5A3418E113F7CF), // -1.800979457880623091361137416183e-6
            F64x2::from_bits(0xBEA25C1238B317FE, 0xBB2C14A72CCA5D1A), // -5.471602987576049206563984667535e-7
            F64x2::from_bits(0x3E57486F6B9C55B0, 0x3ADF4ECAD39D117B), // 2.168393706873286668075385095474e-8
        ],
        k_lo: [
            f64::from_bits(0x3E33E3BF2471A164), // 4.630916600798908e-9
            f64::from_bits(0xBDE78A38A8230039), // -1.7127551119423222e-10
            f64::from_bits(0xBDBE844EA72D8C4E), // -2.775489078897184e-11
            f64::from_bits(0x3D70E24B04393353), // 9.597418489259037e-13
            f64::from_bits(0x3D417979EC729494), // 1.2416390871444032e-13
            f64::from_bits(0xBCF21FDC1C99760D), // -4.024436870032457e-15
            f64::from_bits(0xBCBF08FAE09D7400), // -4.3069822826766145e-16
            f64::from_bits(0x3C6E4EA71C75E7DF), // 1.3143669058584569e-17
            f64::from_bits(0x3C356E99FFE4AF57), // 1.1618330726288413e-18
            f64::from_bits(0xBBE4E76F43496984), // -3.541282373914692e-20
        ],
    },
    // [11.780972450961723, 14.922565104551516)
    ZeroPoly {
        c: f64::from_bits(0x402AA5BAF310E5A2), // 1.3323691936314223e1
        k_hi: [
            F64x2::from_bits(0xBC905DCC62D0D222, 0xB93647158407AF64), // -5.678235636145885426231526782127e-17
            F64x2::from_bits(0x3FCBF3337873A7D9, 0xBC699FEDAB44946D), // 2.183594072478729663959114273106e-1
            F64x2::from_bits(0xBF80C83A2D7ADD23, 0x3C2BA56AD76ABA52), // -8.194403183877490018432881690201e-3
            F64x2::from_bits(0xBFA251858011816B, 0x3C2DE4414AD41E14), // -3.577820957503060397300497209491e-2
            F64x2::from_bits(0x3F559EB160BF72CE, 0xBBF25D2A62A3B155), // 1.319573612810365424255753266496e-3
            F64x2::from_bits(0x3F5C5BCE33AF2D77, 0x3BED6D9FF14FF4D0), // 1.730872506174955568582388547321e-3
            F64x2::from_bits(0xBF10413E306E0037, 0xBB7152DEE66C23AD), // -6.200735161421352569592170276501e-5
            F64x2::from_bits(0xBF04A6704D05ACE7, 0xBBA52996C76931FB), // -3.938703742075740025265611291015e-5
            F64x2::from_bits(0x3EB6C43EEDFED993, 0x3B53CB22E5439A0E), // 1.356994225416156885243449159362e-6
            F64x2::from_bits(0x3EA16ABD7815C6EA, 0xBB3AFF450560C3C8), // 5.190656652677219986986156814340e-7
            F64x2::from_bits(0xBE5257F16F5E8ADD, 0x3AF79A2BFDAD574B), // -1.708374149752036405746231985218e-8
        ],
        k_lo: [
            f64::from_bits(0xBE332DB1B4A162E7), // -4.465340828075029e-9
            f64::from_bits(0x3DE33ACCD029C8C6), // 1.3991439264255694e-10
            f64::from_bits(0x3DBDC8F5580973B1), // 2.7089293917679014e-11
            f64::from_bits(0xBD6C6513A69186D4), // -8.070296387144904e-13
            f64::from_bits(0xBD41358115EA334F), // -1.2227730607616847e-13
            f64::from_bits(0x3CEF325725478310), // 3.46352305021452e-15
            f64::from_bits(0x3CBEC421BD597247), // 4.2696596017275916e-16
            f64::from_bits(0xBC6A925F61E39308), // -1.152366740208341e-17
            f64::from_bits(0xBC3557AF85FE7822), // -1.1569804372777888e-18
            f64::from_bits(0x3BE274AB67D565AF), // 3.126524044051488e-20
        ],
    },
    // [14.922565104551516, 18.06415775814131)
    ZeroPoly {
        c: f64::from_bits(0x4030787B360508C5), // 1.6470630050877634e1
        k_hi: [
            F64x2::from_bits(0xBCB6EB905BA2ABFA, 0x395D7D66F88446BC), // -3.180812762837804809051967814367e-16
            F64x2::from_bits(0xBFC925C6FCA08F54, 0x3C563A5CD9337C45), // -1.964653714686571635656939813763e-1
            F64x2::from_bits(0x3F786DD32E059BC4, 0xBC1E13E7360A0F36), // 5.964112206448160822080956619767e-3
            F64x2::from_bits(0x3FA09463BBD0367E, 0x3C44BD78A9DB1FBD), // 3.238212268489081227574042481533e-2
            F64x2::from_bits(0xBF4FDA0298C87701, 0xBBB77BFB4D2CFCF2), // -9.720337562304227653603720730645e-4
            F64x2::from_bits(0xBF59F4BE60758FAF, 0x3BD149C018AC71DC), // -1.584230341756569700459573776507e-3
            F64x2::from_bits(0x3F0877991AF9D1FD, 0x3B850A48E825F646), // 4.666744223139408197630338630309e-5
            F64x2::from_bits(0x3F032CB00EE8C1CC, 0x3B7CEE35543DBFDF), // 3.657257353884314272870031934485e-5
            F64x2::from_bits(0xBEB19D8CE8C361BB, 0x3B543E9D66EFB1CE), // -1.049961528434586383274735870666e-6
            F64x2::from_bits(0xBEA06A042FBB8C64, 0xBB45DC432C9D4643), // -4.891790860963335797114686929770e-7
            F64x2::from_bits(0x3E4D3A689E6982F4, 0xBAD61505C69CA8A9), // 1.361042207859581925712583342752e-8
        ],
        k_lo: [
            f64::from_bits(0x3E325108C4BC5275), // 4.264651808578063e-9
            f64::from_bits(0xBDDF7B8E9BC9F9E4), // -1.1453319451463557e-10
            f64::from_bits(0xBDBCC40CF520C3A9), // -2.6162363352555305e-11
            f64::from_bits(0x3D67CD773A17FE42), // 6.765105953421137e-13
            f64::from_bits(0x3D40C582D57E47CF), // 1.1916887116405969e-13
            f64::from_bits(0xBCEAB00C1CBC9459), // -2.9629282159222432e-15
            f64::from_bits(0xBCBE33B424A74761), // -4.1913649657590727e-16
            f64::from_bits(0x3C671F81AAA9E9FB), // 1.0028034091542494e-17
            f64::from_bits(0x3C35110E02562AC1), // 1.1420237592141152e-18
            f64::from_bits(0xBBE037BB971CBFC0), // -2.7473862876808616e-20
        ],
    },
    // [18.06415775814131, 21.205750411731103)
    ZeroPoly {
        c: f64::from_bits(0x40339DA8E7416CA4), // 1.9615858510468243e1
        k_hi: [
            F64x2::from_bits(0x3CAA10B2F7B4E69D, 0xB93C423F01E04FCC), // 1.808638713027252511465600506767e-16
            F64x2::from_bits(0x3FC70C511227D5AA, 0xBC3345CD370C6C9D), // 1.800633753443155449890955813952e-1
            F64x2::from_bits(0xBF72CCB0E9755942, 0x3C1C7A7B99F27A2D), // -4.589739858906151342824096953957e-3
            F64x2::from_bits(0xBF9E7DC08E70E999, 0x3C1C48BCC7814ECA), // -2.977658147231352112586614763530e-2
            F64x2::from_bits(0x3F48ACDC5B058C52, 0xBBAEE595764F29AE), // 7.530284838533996048703382873984e-4
            F64x2::from_bits(0x3F580503724AD308, 0x3BFD06DBD9EFE830), // 1.466039052656133602295640284648e-3
            F64x2::from_bits(0xBF032EE4CA1FCB23, 0xBB8C7824BB9D646C), // -3.658900938310050343838867394925e-5
            F64x2::from_bits(0xBF01E5D2836C8D73, 0x3BA7C060C1A1CC5D), // -3.413723646008671967952927070826e-5
            F64x2::from_bits(0x3EAC129F077BB542, 0x3B33B66034D590E8), // 8.366328207700393695954442448666e-7
            F64x2::from_bits(0x3E9EF161591151D4, 0xBB2209DBE178843B), // 4.610850164500622480713673745170e-7
            F64x2::from_bits(0xBE47B9BB07F34E99, 0xBACB8F41F14FD93B), // -1.104805158331828744617389512323e-8
        ],
        k_lo: [
            f64::from_bits(0xBE316F393747926F), // -4.059278125850526e-9
            f64::from_bits(0x3DDA0BC8673EAEA6), // 9.475489300415667e-11
            f64::from_bits(0x3DBBA135E96AD3BF), // 2.5129092026486592e-11
            f64::from_bits(0xBD640D54848D5E05), // -5.699141322852765e-13
            f64::from_bits(0xBD403D00FFD2FF93), // -1.1538003617993764e-13
            f64::from_bits(0x3CE6DB77FE348629), // 2.537670040515241e-15
            f64::from_bits(0x3CBD713898EB1536), // 4.085935736534707e-16
            f64::from_bits(0xBC6418D3117440CE), // -8.7156716934045e-18
            f64::from_bits(0xBC34A6ADC1074787), // -1.1194977657510422e-18
            f64::from_bits(0x3BDC78F154064D89), // 2.4117088311299257e-20
        ],
    },
    // [21.205750411731103, 24.347343065320896)
    ZeroPoly {
        c: f64::from_bits(0x4036C294E3D4D8AC), // 2.2760084380592772e1
        k_hi: [
            F64x2::from_bits(0xBC97BC6D5A660382, 0xB9258F5AF1ED4F1B), // -8.235094410619702873629994596085e-17
            F64x2::from_bits(0xBFC5664E13B70621, 0xBC6D607BC39364A0), // -1.671846004738180418717730832219e-1
            F64x2::from_bits(0x3F6E16555E108E25, 0xBC01E6B9C11D0ED6), // 3.672758801728697557314248604842e-3
            F64x2::from_bits(0x3F9C5E1AD9FB2F40, 0xBC2C2B366BAF731A), // 2.770273166133496646455511532751e-2
            F64x2::from_bits(0xBF43D369F958E589, 0x3BCF6A9AF90EB180), // -6.050364924653921138395885692908e-4
            F64x2::from_bits(0xBF566F4EC27A96E9, 0x3BF465BEF0FD5E88), // -1.369311250451083423267871236665e-3
            F64x2::from_bits(0x3EFF0DE0532652FD, 0xBB95D43FCF6B82DC), // 2.961559693799834919606842060099e-5
            F64x2::from_bits(0x3F00CF2643414079, 0xBB81C0E273A221B0), // 3.206096190011635605156577405333e-5
            F64x2::from_bits(0xBEA6F46D51E579A2, 0xBB46CF105881396E), // -6.841061439935231914535880807328e-7
            F64x2::from_bits(0xBE9D407F7C245E05, 0x3B0EFB763915226E), // -4.358879516371235679168872952933e-7
            F64x2::from_bits(0x3E43A33CD9E0D11E, 0x3ADE90ECD7CA31AF), // 9.144492104601078788141447364585e-9
        ],
        k_lo: [
            f64::from_bits(0x3E309901B09684CF), // 3.864448991197145e-9
            f64::from_bits(0xBDD5D856A6431365), // -7.9471909578743e-11
            f64::from_bits(0xBDBA7CBCC904F961), // -2.4090018663660413e-11
            f64::from_bits(0x3D610B62FEC019A2), // 4.844332379330664e-13
            f64::from_bits(0x3D3F56979C316DE4), // 1.1133583259586171e-13
            f64::from_bits(0xBCE3AE7DD4C84C33), // -2.185097384305494e-15
            f64::from_bits(0xBCBC9460CB5E0F17), // -3.966216517075168e-16
            f64::from_bits(0x3C6183AB9C184C2C), // 7.595633023118941e-18
            f64::from_bits(0x3C342602FF3784DD), // 1.092251464495832e-18
            f64::from_bits(0xBBD909EC6AB7D13D), // -2.1208657799802166e-20
        ],
    },
    // [24.347343065320896, 27.48893571891069)
    ZeroPoly {
        c: f64::from_bits(0x4039E7570DCEA106), // 2.5903672087618382e1
        k_hi: [
            F64x2::from_bits(0xBC961C29FAC28FDF, 0x39300BBF534F6FAD), // -7.670952608123098843713855097538e-17
            F64x2::from_bits(0x3FC40F90793605BB, 0x3C69EB027F4C118C), // 1.567249862528522331974983081318e-1
            F64x2::from_bits(0xBF68C833077FBE56, 0x3BA541E72E7459DF), // -3.025149981105627707525777142841e-3
            F64x2::from_bits(0xBF9AA0CE0421D1A8, 0xBC2EACC20E9F9E41), // -2.600404644222612054893213912486e-2
            F64x2::from_bits(0x3F405FA598EF5CFF, 0x3BEFFE31965509B6), // 4.996832448068976874745546077682e-4
            F64x2::from_bits(0x3F551D30D78AB527, 0xBBF9D231765AA647), // 1.288697907664124218549535447909e-3
            F64x2::from_bits(0xBEF9C5807675C5DD, 0xBB9DFBF9A07334D2), // -2.457760947777674300578170319492e-5
            F64x2::from_bits(0xBEFFC1BBF57E3A9B, 0x3B8EDA0D892493D4), // -3.028561999901067645845902798706e-5
            F64x2::from_bits(0x3EA32DFEA2518F83, 0x3B33FC3FEAC302A8), // 5.715986090112964645563078793056e-7
            F64x2::from_bits(0x3E9BC212085D9D5C, 0x3B363547120D4262), // 4.136277386074759336347839839961e-7
            F64x2::from_bits(0xBE408B946D65FBB5, 0x3AEEC2B87D731C92), // -7.704474763680412448189209693745e-9
        ],
        k_lo: [
            f64::from_bits(0xBE2FA8F9D8B7F3B8), // -3.6857163482240288e-9
            f64::from_bits(0x3DD293FE1552133A), // 6.758671862636771e-11
            f64::from_bits(0x3DB96544BBCC6A0F), // 2.309714549602971e-11
            f64::from_bits(0xBD5D4750DA6BF4FC), // -4.1607361072341933e-13
            f64::from_bits(0xBD3E340F3E9EBD5B), // -1.073038817419678e-13
            f64::from_bits(0x3CE112BC942A5E83), // 1.8955048623698924e-15
            f64::from_bits(0x3CBBAE6180EBA3B1), // 3.8415347690109557e-16
            f64::from_bits(0xBC5EAE41C2D6EED2), // -6.6528143529371084e-18
            f64::from_bits(0xBC339A12DC6330E8), // -1.0626184336918516e-18
            f64::from_bits(0x3BD61A02860948F8), // 1.8720784239796538e-20
        ],
    },
    // [27.48893571891069, 30.630528372500482)
    ZeroPoly {
        c: f64::from_bits(0x403D0BFCF471FCCC), // 2.9046828534916855e1
        k_hi: [
            F64x2::from_bits(0xBC87E3B01386785F, 0x392A842CF25B5D32), // -4.144151247975630458445940249352e-17
            F64x2::from_bits(0xBFC2F2072E638CF4, 0x3C6BE755456303E0), // -1.480111099727775426650330432997e-1
            F64x2::from_bits(0x3F64DF208BBD4521, 0x3BE912224212909E), // 2.547801557661558612277795277440e-3
            F64x2::from_bits(0x3F992BB5E1E159FC, 0xBC34DD401D5F79CA), // 2.458080474056066497960651552610e-2
            F64x2::from_bits(0xBF3BA181C06897EC, 0xBBDE94C5B024C6BF), // -4.216138626439425682524523930784e-4
            F64x2::from_bits(0xBF53FE9D5BAA4A3D, 0x3BE62772A7ED0B23), // -1.220372838962049266984295942964e-3
            F64x2::from_bits(0x3EF5D17602B01CC2, 0xBB8AE18AA3C5F6C1), // 2.080746359569836823069324767543e-5
            F64x2::from_bits(0x3EFE26D3747FE7E2, 0x3B348D44EFDA8352), // 2.875486759987961711557461921876e-5
            F64x2::from_bits(0xBEA0509768AB7128, 0x3B26BEC09F6F1130), // -4.862192367299928933272408718257e-7
            F64x2::from_bits(0xBE9A70F232D9A2FE, 0xBB2D3688C4EDB877), // -3.940045186558046159775092868588e-7
            F64x2::from_bits(0x3E3C509252E04D70, 0xBACF26279C4946A2), // 6.592537444785351141334278501832e-9
        ],
        k_lo: [
            f64::from_bits(0x3E2E454FEDE551A3), // 3.5239792049430794e-9
            f64::from_bits(0xBDD0015B06B8BE87), // -5.822692474047173e-11
            f64::from_bits(0xBDB860E94B8F7E5A), // -2.2172170976418437e-11
            f64::from_bits(0x3D59691EE909A84A), // 3.611067413629974e-13
            f64::from_bits(0x3D3D1CDEF981A18C), // 1.0342936222297717e-13
            f64::from_bits(0xBCDDDCF14F10DE74), // -1.6577326778432272e-15
            f64::from_bits(0xBCBACA1ED43DE40E), // -3.717794532546538e-16
            f64::from_bits(0x3C5B09D1A64D3935), // 5.8630087003376094e-18
            f64::from_bits(0x3C330AAB506FC955), // 1.0322513540189964e-18
            f64::from_bits(0xBBD39AB1DD200A95), // -1.6605468090513646e-20
        ],
    },
    // [30.630528372500482, 33.772121026090275)
    ZeroPoly {
        c: f64::from_bits(0x404018476E6B2BF0), // 3.2189679910974405e1
        k_hi: [
            F64x2::from_bits(0x3CAF5EFD41F756B6, 0xB94BD22B49961F02), // 2.176804008263890746530402653210e-16
            F64x2::from_bits(0x3FC1FF5EEC6A01CD, 0xBC4A6C0DF5A2EEE2), // 1.406057981839822472954959134187e-1
            F64x2::from_bits(0xBF61E438B722C4B0, 0x3C00F94D4D1EFD4D), // -2.184019825187087533890837476393e-3
            F64x2::from_bits(0xBF97ED5FFFC1C773, 0x3C08669892707568), // -2.336645124928033914865016086951e-2
            F64x2::from_bits(0x3F37B7997BABDA6F, 0x3BDDFD6FFEB6E190), // 3.618955339982692840399521762850e-4
            F64x2::from_bits(0x3F53081DEF9612C4, 0xBBE1105BABA9A291), // 1.161603197347081123616604248634e-3
            F64x2::from_bits(0xBEF2C5F5EDAFC545, 0xBB9A1E86498CF76C), // -1.790359861425962824860315797545e-5
            F64x2::from_bits(0xBEFCC11A59E136F3, 0xBB684991E4602928), // -2.742224534584606268285136006847e-5
            F64x2::from_bits(0x3E9C2C3A1B80191B, 0xBB3DA94E542E2913), // 4.198068626063494081425677550406e-7
            F64x2::from_bits(0x3E9946D1DAB790A1, 0x3AECD54C5C3D7829), // 3.766512814680661350159287951422e-7
            F64x2::from_bits(0xBE388DB6194895CA, 0x3AC302790AFDD817), // -5.716821145411937627696311446161e-9
        ],
        k_lo: [
            f64::from_bits(0xBE2D04D33BC48DC3), // -3.3782385484001095e-9
            f64::from_bits(0x3DCBE64387CA000D), // 5.074883653383234e-11
            f64::from_bits(0x3DB77142D1F29720), // 2.1320762101760325e-11
            f64::from_bits(0xBD564584C374BE8A), // -3.1649786166445367e-13
            f64::from_bits(0xBD3C15E1013B1566), // -9.977961408544963e-14
            f64::from_bits(0x3CDA5479F668D687), // 1.4616078347977999e-15
            f64::from_bits(0x3CB9EDCAB9AB39B9), // 3.5983541968482367e-16
            f64::from_bits(0xBC57FC2C0F3D6D6B), // -5.200928076900961e-18
            f64::from_bits(0xBC327CD13F15CCC0), // -1.002213061188992e-18
            f64::from_bits(0x3BD17BCA1BD679F5), // 1.480914515895722e-20
        ],
    },
    // [33.772121026090275, 36.91371367968007)
    ZeroPoly {
        c: f64::from_bits(0x4041AA890DC5E97C), // 3.533230755008387e1
        k_hi: [
            F64x2::from_bits(0xBCBF89DCEDB3EA9B, 0x39469B3746B0E6BF), // -4.376849903590193489184104231832e-16
            F64x2::from_bits(0xBFC12DD57BF18AD9, 0xBC43F7436F827837), // -1.342112403100006756126919108679e-1
            F64x2::from_bits(0x3F5F1E1E7F394272, 0xBBF91F7055F5D02E), // 1.899270803637300929440939575747e-3
            F64x2::from_bits(0x3F96D9AFE88301F9, 0x3C1AA220D9D26A76), // 2.231478554381458417900507097818e-2
            F64x2::from_bits(0xBF34A538A48298E9, 0xBBD59C2E3E776FD4), // -3.150237357750900118350624268296e-4
            F64x2::from_bits(0xBF52316250B4AE35, 0xBBE18785B8CF5453), // -1.110406886094217758021127158851e-3
            F64x2::from_bits(0x3EF05F11577B46DC, 0xBB828C108DF92B61), // 1.561294399616315158820874958752e-5
            F64x2::from_bits(0x3EFB86BAD42FC1DA, 0x3B36B9E1086EFE6B), // 2.625111416060415708854700767578e-5
            F64x2::from_bits(0xBE98A1B3A9E92BAE, 0x3B30DCC370867E3B), // -3.670401528584423664312770472042e-7
            F64x2::from_bits(0xBE983DCAF3F8D179, 0xBB15F573EACD6FCD), // -3.612246821411225143350182001209e-7
            F64x2::from_bits(0x3E3589A7CA617E86, 0xBADE9027A5E10E81), // 5.014640402656822759377492962115e-9
        ],
        k_lo: [
            f64::from_bits(0x3E2BE3EE32789893), // 3.2468644611197007e-9
            f64::from_bits(0xBDC8913F1DEAA5E1), // -4.4687784496780064e-11
            f64::from_bits(0xBDB695C377D3CB14), // -2.0540950441620116e-11
            f64::from_bits(0x3D53B25D7A5FC7D9), // 2.7990749432334136e-13
            f64::from_bits(0x3D3B20BBFAE47718), // 9.637754895631752e-14
            f64::from_bits(0xBCD7651D7BEA893C), // -1.2986823361627766e-15
            f64::from_bits(0xBCB91C85207804C6), // -3.484907689399181e-16
            f64::from_bits(0x3C556A5D6F052389), // 4.643743762936763e-18
            f64::from_bits(0x3C31F37D7A19AE45), // 9.731328902623438e-19
            f64::from_bits(0xBBCF5BBDDB761CB9), // -1.3280784477783515e-20
        ],
    },
    // [36.91371367968007, 40.05530633326986)
    ZeroPoly {
        c: f64::from_bits(0x40433CC523D5CB69), // 3.8474766234771614e1
        k_hi: [
            F64x2::from_bits(0xBC9AAAF726A29E97, 0xB932FEB3F4319D7B), // -9.252263442802770616132542889925e-17
            F64x2::from_bits(0x3FC076826CC2C191, 0x3C672CBCDCFBBF39), // 1.286166220720699571751659377652e-1
            F64x2::from_bits(0xBF5B62885E006FF1, 0xBBE263A80285280D), // -1.671441241348322724183181975150e-3
            F64x2::from_bits(0xBF95E7F53001E4B1, 0x3C1CC7AF11D46E58), // -2.139266114725008357942021707116e-2
            F64x2::from_bits(0x3F322EBEB8DC21BB, 0xBBBDA3A5AB0FB034), // 2.774444227360524375460141745739e-4
            F64x2::from_bits(0x3F517444A7A04CD0, 0xBBF2978EDF864984), // 1.065318151181316528599170083272e-3
            F64x2::from_bits(0xBEECE06F1F1FCD3B, 0x3B8D803F57C5CB05), // -1.376948145901409620128498945302e-5
            F64x2::from_bits(0xBEFA7006E6AD9CBD, 0x3B80A1312AF6009E), // -2.521286516401760756183564919792e-5
            F64x2::from_bits(0x3E95C42F02CF18D9, 0x3B24456FFBF9B252), // 3.243437756674553403479314991772e-7
            F64x2::from_bits(0x3E9750CA5E133C7F, 0x3B35AD749B252BA2), // 3.474293333695238559876709652812e-7
            F64x2::from_bits(0xBE3314982DF95EC5, 0x3AD00C3CB5FD27BF), // -4.442512773956995610641574605615e-9
        ],
        k_lo: [
            f64::from_bits(0xBE2ADED752E75E90), // -3.128134771378151e-9
            f64::from_bits(0x3DC5D478489C0BF1), // 3.970846661569224e-11
            f64::from_bits(0x3DB5CCF43EF45FE0), // 1.9827531920706952e-11
            f64::from_bits(0xBD518FCE7AB8D895), // -2.495673978636341e-13
            f64::from_bits(0xBD3A3D63CCCCD104), // -9.322251087394403e-14
            f64::from_bits(0x3CD4EFD314FFCB19), // 1.1622266819023349e-15
            f64::from_bits(0x3CB857915A07A04B), // 3.378139662530934e-16
            f64::from_bits(0xBC533DA5505E8F2A), // -4.172184243181501e-18
            f64::from_bits(0xBC3170488AA93DC6), // -9.453487741780775e-19
            f64::from_bits(0x3BCC46AF6390F864), // 1.1975399971717874e-20
        ],
    },
    // [40.05530633326986, 43.19689898685966)
    ZeroPoly {
        c: f64::from_bits(0x4044CEFCF1734B62), // 4.161709421281445e1
        k_hi: [
            F64x2::from_bits(0x3C9451B6225ACBFB, 0x3933C613C6379F68), // 7.049633584286915508535588689674e-17
            F64x2::from_bits(0xBFBFA8B41711C83A, 0x3C55ACFAD359589B), // -1.236679607698371272289254947113e-1
            F64x2::from_bits(0x3F5857D39699972F, 0x3BD2FCA603DA8EBC), // 1.485783223324554030068524290720e-3
            F64x2::from_bits(0x3F9511C6DADAAA12, 0xBC2E1683C66CDEDE), // 2.057562552242059918431554993903e-2
            F64x2::from_bits(0xBF302C289DBDBD19, 0x3BC273DC599C7BB0), // -2.467726860606186552616047725566e-4
            F64x2::from_bits(0xBF50CC2238D229F9, 0xBBDD636AE71B3CC8), // -1.025231761817678338959449351383e-3
            F64x2::from_bits(0x3EE9B64D5C63665D, 0x3B42D7ACAD809C2A), // 1.226049324680731987237765328690e-5
            F64x2::from_bits(0x3EF976FB023F0F39, 0x3B79EF406FDE73DE), // 2.428509482378592766389453734662e-5
            F64x2::from_bits(0xBE93693BA0B5BD36, 0xBB3188461421BC97), // -2.892474248649520483950282698248e-7
            F64x2::from_bits(0xBE967B9529870BE4, 0x3B3F27810569ECCA), // -3.350190041184976386325738333156e-7
            F64x2::from_bits(0x3E310CB79A2C2C51, 0x3AAFC759AAF30D89), // 3.969687164714230920216444909741e-9
        ],
        k_lo: [
            f64::from_bits(0x3E29F2079F6FAC91), // 3.0204454453968334e-9
            f64::from_bits(0xBDC38D957F59B7B1), // -3.556681331292419e-11
            f64::from_bits(0xBDB51511DB62E968), // -1.91742435398382e-11
            f64::from_bits(0x3D4F8BB546E36517), // 2.241459271778719e-13
            f64::from_bits(0x3D396AF6B7DB7C12), // 9.030226209767321e-14
            f64::from_bits(0xBCD2DC4FC009D22A), // -1.0469731687564795e-15
            f64::from_bits(0xBCB79F20F0969831), // -3.278155021245144e-16
            f64::from_bits(0x3C51636D93E9116B), // 3.770506210673362e-18
            f64::from_bits(0x3C30F3EB14F8E4D9), // 9.190134447430207e-19
            f64::from_bits(0xBBC9A1692D382863), // -1.085494368601435e-20
        ],
    },
    // [43.19689898685966, 46.33849164044945)
    ZeroPoly {
        c: f64::from_bits(0x404661315D6B133F), // 4.475931899765282e1
        k_hi: [
            F64x2::from_bits(0xBCB40032091A4E00, 0xB953CC4AF3EFD26D), // -2.775663515977503714613668548237e-16
            F64x2::from_bits(0x3FBE8727DAA3DAED, 0x3C5D1A2F9814E55E), // 1.192498120106894800951832981801e-1
            F64x2::from_bits(0xBF55D353E28547B8, 0xBBE1E2A55C764682), // -1.332122725291358072742094457347e-3
            F64x2::from_bits(0xBF94524D4813CC25, 0xBC329396E7DEF0FA), // -1.984520676064617964679726962932e-2
            F64x2::from_bits(0x3F2D037574E281C7, 0x3BA9B67F41A9F545), // 2.213555221191166997869545574352e-4
            F64x2::from_bits(0x3F50356BB747A764, 0xBBDC2A876D040EF5), // 9.892990029762510737054148559188e-4
            F64x2::from_bits(0xBEE7156BFCCEF29B, 0x3B8ADFE268538DF4), // -1.100715589951702613246287475930e-5
            F64x2::from_bits(0xBEF896D7DC819F72, 0xBB9CC20F22F4B789), // -2.345011833462758768140798654286e-5
            F64x2::from_bits(0x3E9172C6DADF4369, 0x3B379FB5C554B854), // 2.600006280247789630489395673775e-7
            F64x2::from_bits(0x3E95BAAE8EFC060B, 0x3B3BDEA3B523018D), // 3.237907000159509912339348534232e-7
            F64x2::from_bits(0xBE2EB347EB4FC67B, 0x3ABAFCC1A97CE8BB), // -3.573987184754700487379127660494e-9
        ],
        k_lo: [
            f64::from_bits(0xBE291A60A70C5499), // -2.9223781661825006e-9
            f64::from_bits(0x3DC1A345AA455722), // 3.208293811508972e-11
            f64::from_bits(0x3DB46C56A2981175), // 1.857478941675978e-11
            f64::from_bits(0xBD4C84BB9A08AFAE), // -2.0263604181062368e-13
            f64::from_bits(0xBD38A836D5C0B70A), // -8.759956925270511e-14
            f64::from_bits(0x3CD117A8566E57B1), // 9.488194884362298e-16
            f64::from_bits(0x3CB6F2D07B9D692E), // 3.1847432602408713e-16
            f64::from_bits(0xBC4F9984AE252846), // -3.4260442562383595e-18
            f64::from_bits(0xBC307E954EEC6F8F), // -8.941667807540163e-19
            f64::from_bits(0x3BC758E499E2CD35), // 9.887939981327252e-21
        ],
    },
    // [46.33849164044945, 49.480084294039244)
    ZeroPoly {
        c: f64::from_bits(0x4047F36312028AD6), // 4.790146088718545e1
        k_hi: [
            F64x2::from_bits(0xBCBCCB5A05A6E4AA, 0xB95FFDFC496E916E), // -3.996017735860100430037733577491e-16
            F64x2::from_bits(0xBFBD8293AA55D18F, 0x3C438D934F192149), // -1.152736941201679527929929836829e-1
            F64x2::from_bits(0x3F53B6BEB83F292E, 0x3BE303FC8ECF5BD2), // 1.203237771721331631924979325028e-3
            F64x2::from_bits(0x3F93A5CCBC12A67B, 0xBC3C77973D90076E), // 1.918716333299251062816367189493e-2
            F64x2::from_bits(0xBF2A3765D26AA68F, 0x3BCCD7E7EB418B42), // -2.000152391777209939328188653821e-4
            F64x2::from_bits(0xBF4F5AB33748C213, 0x3BE76BA4131883DE), // -9.568572118663929645742630667992e-4
            F64x2::from_bits(0x3EE4DF6F1C257BA6, 0x3B8ACE63EDD6658E), // 9.952921463605491692048612616094e-6
            F64x2::from_bits(0x3EF7CBD49C315BA1, 0xBB7C451EF899EA0B), // 2.269383709249293528787195269356e-5
            F64x2::from_bits(0xBE8F96098CF077AF, 0x3B2DD78C9519B700), // -2.353346588443222304393188835494e-7
            F64x2::from_bits(0xBE950B37DD432BED, 0xBB242D1856B4D65C), // -3.135773714134773136892353251204e-7
            F64x2::from_bits(0x3E2BD2E6405E88FF, 0x3AC0AA31E11429D5), // 3.239119641738049601551203567941e-9
        ],
        k_lo: [
            f64::from_bits(0x3E285530DEF03749), // 2.8327080607455895e-9
            f64::from_bits(0xBDC0029E22237D49), // -2.9122430393713905e-11
            f64::from_bits(0xBDB3D11ADE7A0DC0), // -1.8023289375527938e-11
            f64::from_bits(0x3D49EF07D7C0041D), // 1.8427011682738686e-13
            f64::from_bits(0x3D37F3CAA8AE3D17), // 8.509570321116314e-14
            f64::from_bits(0xBCCF2637AAD6D1B7), // -8.645663883756355e-16
            f64::from_bits(0xBCB651F1BB918441), // -3.097535394585205e-16
            f64::from_bits(0x3C4CDB5B4D209EF0), // 3.1286672359795462e-18
            f64::from_bits(0x3C301027ACFFFDDF), // 8.707826888143537e-19
            f64::from_bits(0xBBC55D908C1BE782), // -9.048635655142573e-21
        ],
    },
];

// GENERATE: bessel_table y0 16 3 0.5 11 10
// Y0_TBL[i] has the coefficients `K0, ..., K20` such that
// f(c + d) ~= K0 + K1 * d + ... + K20 * d^20, where `c` is very close to the
// (i + 1)-th positive zero of y0, for the values of the interval of that
// zero, [(i - 0.25) * π, (i + 0.75) * π) (the first one starts at 0.5)
// f(x) = y0(x) - 2/π * j0(x) * ln(x / c) in the first 3 intervals
// f(x) = y0(x) in the other ones
static Y0_TBL: [ZeroPoly; 16] = [
    // [0.5, 2.356194490192345)
    ZeroPoly {
        c: f64::from_bits(0x3FEC982EB8D417EA), // 8.935769662791675e-1
        k_hi: [
            F64x2::from_bits(0xBC7AF74BFA0F1304, 0xB9111721BF3251BD), // -2.338927928406210311869215393415e-17
            F64x2::from_bits(0x3FD3582AFCD83AFC, 0x3C6C8C30C52D1E40), // 3.022563428102527218215442195571e-1
            F64x2::from_bits(0x3FBE52DB13E9C881, 0xBC5F7E99692B732F), // 1.184517787567376316245903223178e-1
            F64x2::from_bits(0xBFAB764E92760B01, 0x3C24E45573A9E05C), // -5.363698518615046066931111557545e-2
            F64x2::from_bits(0xBF86205E65D73516, 0x3C293A36B36EB0F7), // -1.080392597477923128834381570756e-2
            F64x2::from_bits(0x3F65E13DDA4F97AC, 0x3C0E6EC982997C17), // 2.670880132870720317369395386494e-3
            F64x2::from_bits(0x3F37A37794B3566A, 0xBBBFFEC414DC9E42), // 3.606955476248343225072521477063e-4
            F64x2::from_bits(0xBF105CEA95501129, 0xBB943A36F52F0120), // -6.241971762227918768850420427512e-5
            F64x2::from_bits(0xBEDA8E17A0608B80, 0xBB56CBB285083B61), // -6.331216814844660325174950288441e-6
            F64x2::from_bits(0x3EAC72852F4E2AA0, 0xBB4588D4CA9DB932), // 8.477969389725677122062313791373e-7
            F64x2::from_bits(0x3E7278525D2C875E, 0xBB1019FF4AD5859F), // 6.880613704883587463985334955830e-8
        ],
        k_lo: [
            f64::from_bits(0xBE402D626FE033E3), // -7.533134557141433e-9
            f64::from_bits(0xBE018065614504E4), // -5.093620547389035e-10
            f64::from_bits(0x3DC9F66434FAD370), // 4.722545152025803e-11
            f64::from_bits(0x3D8810E57ED1C79E), // 2.735987643714419e-12
            f64::from_bits(0xBD4EFB4DF42BD1D0), // -2.2013792197479403e-13
            f64::from_bits(0xBD091D40696D86BF), // -1.1152973698575505e-14
            f64::from_bits(0x3CCC83DD9E426C44), // 7.914530246561629e-16
            f64::from_bits(0x3C8521561EC527CA), // 3.665508927619157e-17
            f64::from_bits(0xBC481D1BB2F93C44), // -2.6144130160431015e-18
            f64::from_bits(0xBBDC8945ADFC6BE1), // -2.4171118086946664e-20
        ],
    },
    // [2.356194490192345, 5.497787143782138)
    ZeroPoly {
        c: f64::from_bits(0x400FA9534D98569C), // 3.957678419314858e0
        k_hi: [
            F64x2::from_bits(0xBC88FA8956B4B481, 0xB9254CFBECD28BBA), // -4.333106464293519638851066383110e-17
            F64x2::from_bits(0xBFD5A61DC95DB6AB, 0xBC65F81C4DD43EE5), // -3.382639376130593495758521459027e-1
            F64x2::from_bits(0x3FA1C77746750248, 0x3C34BE5541DE9637), // 3.472492919036490711988745612607e-2
            F64x2::from_bits(0x3F9E7DB187E2ACF3, 0x3C1258BCEC0C8DD6), // 2.977635757330872001586765548317e-2
            F64x2::from_bits(0x3F4118CD1943FBDD, 0x3BEE1FE2BBFCEE36), // 5.217553575959062195340439419102e-4
            F64x2::from_bits(0xBF5402F6BD729ACF, 0xBBEF0BBEA119B0AB), // -1.221409756718184166434936918337e-3
            F64x2::from_bits(0xBF1545E4BB2E3A41, 0xBBA901D4A243A7CE), // -8.115013660712438552591207254089e-5
            F64x2::from_bits(0x3EFCE055CB0E4D01, 0x3B9D405BF6113F90), // 2.753859434141432493589430976930e-5
            F64x2::from_bits(0x3EC19F9830EEFAA4, 0xBB30E58286F309FD), // 2.100874901128065420581620576641e-6
            F64x2::from_bits(0xBE99610D7665F252, 0x3B3F3C048D39C52D), // -3.781782339713302003089386693693e-7
            F64x2::from_bits(0xBE5DCAA92D24638C, 0xBAE1687234EC32D7), // -2.774563052322065308939371267056e-8
        ],
        k_lo: [
            f64::from_bits(0x3E2D906A98D6382C), // 3.4417173058520393e-9
            f64::from_bits(0x3DEFCFB749F586A1), // 2.3145832904777095e-10
            f64::from_bits(0xBDB859FCC09FF07B), // -2.2147572004700428e-11
            f64::from_bits(0xBD77BA43D5F700B8), // -1.3487577684922114e-12
            f64::from_bits(0x3D3DCAE82327F83F), // 1.0584459401219434e-13
            f64::from_bits(0x3CFA4839C98FE4B9), // 5.83580556434992e-15
            f64::from_bits(0xBCBC18EF27A8E86B), // -3.8992974425020053e-16
            f64::from_bits(0xBC768B9EF2B8133E), // -1.955501221528725e-17
            f64::from_bits(0x3C3474C3314FFCCD), // 1.1089275875124503e-18
            f64::from_bits(0x3BEEAEB443DF6754), // 5.197807199226576e-20
        ],
    },
    // [5.497787143782138, 8.63937979737193)
    ZeroPoly {
        c: f64::from_bits(0x401C581DC4E72103), // 7.086051060301773e0
        k_hi: [
            F64x2::from_bits(0x3C7E91B198D39CE2, 0x391EE1B6E264023D), // 2.651448172550504090177680043808e-17
            F64x2::from_bits(0x3FD17C223FD00A15, 0xBC7B85E15EB8C82E), // 2.732015250368544961841794539380e-1
            F64x2::from_bits(0xBF91CE4B4929E6E3, 0x3C304BD6AEE0B03A), // -1.738851198249492366807992679075e-2
            F64x2::from_bits(0xBFA1EE4A83DD5E20, 0xBBFFD590185D09FD), // -3.502114162136771302286137875676e-2
            F64x2::from_bits(0x3F582C42A99BCDE7, 0xBBF9DB157EB6CE69), // 1.475396251802924039477760485535e-3
            F64x2::from_bits(0x3F5533A61F5E373A, 0xBBC6798143A584F3), // 1.294052342447794900550228134605e-3
            F64x2::from_bits(0xBEEC496B3AE0374D, 0xBB664C94590DA613), // -1.348819372639820324332072233021e-5
            F64x2::from_bits(0xBEFA9BACEA182936, 0x3B9656CFB54F0ECC), // -2.537546845896820333322613354062e-5
            F64x2::from_bits(0xBE982184CD1DC0D0, 0xBB26F6B5E052A7AA), // -3.595789169674497873425036898839e-7
            F64x2::from_bits(0x3E95294E3EB1F5A1, 0xBB37948DF49222DD), // 3.153286899994205218675915170943e-7
            F64x2::from_bits(0x3E43D47AC2A6EDE2, 0x3AEDC25CD24CF898), // 9.234062477216568347953378297686e-9
        ],
        k_lo: [
            f64::from_bits(0xBE272D0C30D0818D), // -2.69803768780994e-9
            f64::from_bits(0xBDDB81399FA295A7), // -1.0006182680070107e-10
            f64::from_bits(0x3DB26B7E98189B60), // 1.6752801841955378e-11
            f64::from_bits(0x3D67C1F10C7CF402), // 6.752311595196467e-13
            f64::from_bits(0xBD361924BCAD627A), // -7.850863714668582e-14
            f64::from_bits(0xBCECF1C28263F0B8), // -3.213471070012809e-15
            f64::from_bits(0x3CB4A2E67192A1BA), // 2.863865920930908e-16
            f64::from_bits(0x3C6A8C34E023FD66), // 1.1513221720113605e-17
            f64::from_bits(0xBC2DEC640539A0C9), // -8.110754144682401e-19
            f64::from_bits(0xBBE2D557345A9E96), // -3.1904956048636816e-20
        ],
    },
    // [8.63937979737193, 11.780972450961723)
    ZeroPoly {
        c: f64::from_bits(0x402471D735A47D58), // 1.0222345043496418e1
        k_hi: [
            F64x2::from_bits(0xBCACABD7C823BF21, 0x394F241D0B14C882), // -1.989468403830483713210225563876e-16
            F64x2::from_bits(0xBFCFF635CC72B9F0, 0xBC236B301A6E8E76), // -2.497012375146847601149351609156e-1
            F64x2::from_bits(0x3F89036451FF57FD, 0x3C279756891D526C), // 1.221350074039761578299074418454e-2
            F64x2::from_bits(0x3FA4E667A71556AE, 0x3C4C0CE8B1BAAA02), // 4.082034983245584918087820256260e-2
            F64x2::from_bits(0xBF60325EE41E911E, 0x3C0EA83AE9CEF368), // -1.977143606341282313236376467946e-3
            F64x2::from_bits(0xBF5FE23914FB9125, 0xBBF1CB199B711B0B), // -1.946025604348264082936208272724e-3
            F64x2::from_bits(0x3F17F84D7C50E4D4, 0x3B8FC5B4ACC1A891), // 9.143803534152349215894981410704e-5
            F64x2::from_bits(0x3F06AFDD57BE1DEB, 0x3B8FD460F90A691B), // 4.327196344805133996774680593823e-5
            F64x2::from_bits(0xBEC04053ABF437D9, 0x3B647107939441F1), // -1.937303153127183098325069423028e-6
            F64x2::from_bits(0xBEA2AEA9EC48C02A, 0x3B491D6E7F14B018), // -5.567753410037251304930531276099e-7
            F64x2::from_bits(0x3E593EB9F1DD1DC9, 0xBAF97991FC27D6DF), // 2.351126148568337191297350391927e-8
        ],
        k_lo: [
            f64::from_bits(0x3E3428A3A3D0A0FC), // 4.693574026867159e-9
            f64::from_bits(0xBDE99D8C42A0555B), // -1.8637672075205656e-10
            f64::from_bits(0xBDBEC9803F97692F), // -2.800071630676549e-11
            f64::from_bits(0x3D7261363AE498BE), // 1.044766903243219e-12
            f64::from_bits(0x3D418F2B7B240DED), // 1.247660266137669e-13
            f64::from_bits(0xBCF3BC8F5FC0CD12), // -4.382397270445356e-15
            f64::from_bits(0xBCBF0A9182D1C824), // -4.3078433628415006e-16
            f64::from_bits(0x3C706ADD712BCA00), // 1.4239860538902475e-17
            f64::from_bits(0x3C35662560AFB500), // 1.160042539357893e-18
            f64::from_bits(0xBBE52D783CB51500), // -3.58762768623728e-20
        ],
    },
    // [11.780972450961723, 14.922565104551516)
    ZeroPoly {
        c: f64::from_bits(0x402AB8E1C4A1E74A), // 1.3361097473872764e1
        k_hi: [
            F64x2::from_bits(0x3CA4D9FE48A91CEF, 0xB9124A9EAD37B762), // 1.446865982116726264470134133718e-16
            F64x2::from_bits(0x3FCBF32A27594007, 0xBC64A9ECBF385B1A), // 2.183582965976713274083468785805e-1
            F64x2::from_bits(0xBF80BC2D84E6523D, 0xBC1698273911428C), // -8.171420686985725976636450055039e-3
            F64x2::from_bits(0xBFA26CAB38A8B367, 0xBC21B4FB746A6649), // -3.598532740247755621598599250358e-2
            F64x2::from_bits(0x3F55F03E47165D8D, 0x3BE2233DF5FE1EA2), // 1.339016740848328409206201793415e-3
            F64x2::from_bits(0x3F5CAAA76E34992C, 0xBBFD2EBBB25C979A), // 1.749671465242362283434208065293e-3
            F64x2::from_bits(0xBF10C5F18C46D03D, 0x3B5A84EF9F7CF5C9), // -6.398474494667765694030662521220e-5
            F64x2::from_bits(0xBF04F0AF7D46CD20, 0x3B9399F81459AFCE), // -3.994021941213844115495848701475e-5
            F64x2::from_bits(0x3EB7D1E28094E12A, 0xBB16A14249B7BFCD), // 1.419774436973713454850484399927e-6
            F64x2::from_bits(0x3EA1AD0731226B38, 0x3B3FABC32E12462D), // 5.267826017786499602339646726068e-7
            F64x2::from_bits(0xBE536021C8CABFA2, 0x3AF1B1FE0C45E98E), // -1.804485498757483031747466437464e-8
        ],
        k_lo: [
            f64::from_bits(0xBE3371AE6C62D7B0), // -4.527174807063323e-9
            f64::from_bits(0x3DE46CE07760956C), // 1.4861361766954628e-10
            f64::from_bits(0x3DBE235D246E7568), // 2.7410478641166424e-11
            f64::from_bits(0xBD6E408DF240C9B5), // -8.598182696649264e-13
            f64::from_bits(0xBD415F02E80369EE), // -1.234293598646058e-13
            f64::from_bits(0x3CF0A4443C19B3DB), // 3.6951921922144316e-15
            f64::from_bits(0x3CBEFA4B191523B7), // 4.299020640060843e-16
            f64::from_bits(0xBC6C4F1F4E9C8D6A), // -1.2277102710356238e-17
            f64::from_bits(0xBC356DCC95B55B35), // -1.161663157501424e-18
            f64::from_bits(0x3BE273662873C163), // 3.1256832987531407e-20
        ],
    },
    // [14.922565104551516, 18.06415775814131)
    ZeroPoly {
        c: f64::from_bits(0x4030803C74003214), // 1.650092244152809e1
        k_hi: [
            F64x2::from_bits(0x3CACD827F6C074A2, 0xB90ECAB22F071784), // 2.001479523805861127488107404913e-16
            F64x2::from_bits(0xBFC925C35988EE2A, 0x3C5C9C8BDB78B931), // -1.964649378950167828559332661674e-1
            F64x2::from_bits(0x3F786254936758AB, 0x3C1C3999B540D17C), // 5.953150152398982639982391197984e-3
            F64x2::from_bits(0x3FA0A4512039D6A3, 0xBC013EED97025C8D), // 3.250363844283057206923083437857e-2
            F64x2::from_bits(0xBF5013B38CFB926C, 0xBBF3FB300D8059FA), // -9.812596726358714990554479698938e-4
            F64x2::from_bits(0xBF5A24A1215F6685, 0xBBF38387CA36A756), // -1.595647195691192217330404524060e-3
            F64x2::from_bits(0x3F08F91421377F84, 0xBB9867D7A5AC9413), // 4.763214760905492228849752590903e-5
            F64x2::from_bits(0x3F035D17CEC0170A, 0xBBA3A6D40565B02F), // 3.693322091184558415680372188825e-5
            F64x2::from_bits(0xBEB2283A93113FAC, 0x3B5B0A40D8F686A0), // -1.082250104884173611192565637152e-6
            F64x2::from_bits(0xBEA099E71392DC28, 0x3B0BA53C4B509613), // -4.947537838505266026529553658040e-7
            F64x2::from_bits(0x3E4E5DE01E2D94DC, 0xBAE1BD286F79D24F), // 1.414059708709926719512229300024e-8
        ],
        k_lo: [
            f64::from_bits(0x3E32885854D7E361), // 4.3149566814520046e-9
            f64::from_bits(0xBDE0730C294987A9), // -1.19685168584541e-10
            f64::from_bits(0xBDBD1743BACD83BE), // -2.645799870441007e-11
            f64::from_bits(0x3D68F8725C1DB843), // 7.097041530101967e-13
            f64::from_bits(0x3D40F1620B9078AE), // 1.2038656152396206e-13
            f64::from_bits(0xBCEC10D38419A05C), // -3.1159216849232407e-15
            f64::from_bits(0xBCBE77AFFFD82F70), // -4.228219064817534e-16
            f64::from_bits(0x3C684F4F73992042), // 1.0542697827816257e-17
            f64::from_bits(0x3C3536CEB7590C1B), // 1.1500182175327244e-18
            f64::from_bits(0xBBE019A3A673AE21), // -2.727472069323771e-20
        ],
    },
    // [18.06415775814131, 21.205750411731103)
    ZeroPoly {
        c: f64::from_bits(0x4033A42CDF5FEBD7), // 1.964130970088794e1
        k_hi: [
            F64x2::from_bits(0x3CB1D336C8B0A153, 0x3943A0F5EDE1AF5B), // 2.473723266042011568776349138305e-16
            F64x2::from_bits(0x3FC70C4F66CAB47F, 0xBC501718D9A777A5), // 1.800631763375441750178777380058e-1
            F64x2::from_bits(0xBF72C6731071E9C4, 0xBC151F65392934AD), // -4.583787412338617287371429205877e-3
            F64x2::from_bits(0xBF9E924B85A17360, 0x3C3CAA3C838D4D3C), // -2.985494616619932384298517701550e-2
            F64x2::from_bits(0x3F48D6C364D920DF, 0x3BC33FA593F8F856), // 7.580236503566226521205927547873e-4
            F64x2::from_bits(0x3F58247B02D6B0F5, 0xBBF3CF2DC05FA77D), // 1.473541381851237498500738265997e-3
            F64x2::from_bits(0xBF0376125938564D, 0x3BA6AB1058D9A266), // -3.711932654596999149198374899200e-5
            F64x2::from_bits(0xBF0206DA232A2B42, 0x3BAD2CCFA6372CAB), // -3.438332750701221524023663626839e-5
            F64x2::from_bits(0x3EACAF84DB8713DD, 0xBB4A4A7A43B0D78F), // 8.548981247849362132977333454633e-7
            F64x2::from_bits(0x3E9F36711772D520, 0xBB32906D43C4C305), // 4.651049247504386415827657651462e-7
            F64x2::from_bits(0xBE48654B000F559E, 0xBAEAE89BA1C8B991), // -1.136012173385825801033001953477e-8
        ],
        k_lo: [
            f64::from_bits(0xBE3199EC889F1DEE), // -4.098113967920126e-9
            f64::from_bits(0x3DDAEB1A842C041A), // 9.79284717926094e-11
            f64::from_bits(0x3DBBE67C08CF79F9), // 2.5375202417297064e-11
            f64::from_bits(0xBD64CDF9E6FF438B), // -5.913021384736255e-13
            f64::from_bits(0xBD40648EBFFA9959), // -1.1647787226025597e-13
            f64::from_bits(0x3CE7C82B12D6836F), // 2.6403220999292142e-15
            f64::from_bits(0x3CBDB409B3E5EA05), // 4.122157206113424e-16
            f64::from_bits(0xBC64EA8D113D58F2), // -9.070962304564217e-18
            f64::from_bits(0xBC34D0098C2EE0AE), // -1.1282557833243177e-18
            f64::from_bits(0x3BDC16B9CD2B656D), // 2.3792115844396704e-20
        ],
    },
    // [21.205750411731103, 24.347343065320896)
    ZeroPoly {
        c: f64::from_bits(0x4036C832FD77AC07), // 2.2782028047291558e1
        k_hi: [
            F64x2::from_bits(0x3CB3296749AAF800, 0xB95BDA7FB68ADDA7), // 2.659224548370762811829264430662e-16
            F64x2::from_bits(0xBFC5664D37C37D7B, 0xBC55F24AF4D9E9A1), // -1.671844980510107428877311478269e-1
            F64x2::from_bits(0x3F6E0EE8EC846468, 0x3BF7DB6428B24D18), // 3.669218949778275338020369238582e-3
            F64x2::from_bits(0x3F9C6C415C971B4C, 0xBC30F158FC514BA3), // 2.775671126962757270541692744461e-2
            F64x2::from_bits(0xBF43EC49CB941EF0, 0xBBD8ED5B6318937D), // -6.080017402902687179343808059133e-4
            F64x2::from_bits(0xBF56853B2D047886, 0x3BF8C50DD615AFFC), // -1.374538220324990174563134719846e-3
            F64x2::from_bits(0x3EFF632A9396FCE1, 0xBB5B385BB7C7986D), // 2.993332711808967157505794776581e-5
            F64x2::from_bits(0x3F00E6AFA01AF53C, 0x3B9952FBD30FD7F4), // 3.223632303448817894962731516361e-5
            F64x2::from_bits(0xBEA7541247ACDDE6, 0x3B4A0267F8272665), // -6.952406147236105433293048764239e-7
            F64x2::from_bits(0xBE9D735E8BEB5F1B, 0x3B33B45C91F8EE40), // -4.388490530255317220698659415910e-7
            F64x2::from_bits(0x3E440E6A7CAF27A5, 0x3AE8B7DBE0081F17), // 9.339448233786037960614311452921e-9
        ],
        k_lo: [
            f64::from_bits(0x3E30B9CFAE18A42D), // 3.894284646092304e-9
            f64::from_bits(0xBDD668087E0E5379), // -8.151393406963794e-11
            f64::from_bits(0xBDBAB498D69BD1B2), // -2.428847176598624e-11
            f64::from_bits(0x3D618B5654B835F2), // 4.986386004344308e-13
            f64::from_bits(0x3D3F99D0D1354DCE), // 1.1226874558156846e-13
            f64::from_bits(0xBCE450AA069A798C), // -2.255428553669665e-15
            f64::from_bits(0xBCBCD0660BAF13B9), // -3.998753702175002e-16
            f64::from_bits(0x3C6216600D0DFD7E), // 7.844160703740035e-18
            f64::from_bits(0x3C344D92D340D227), // 1.1006290076855078e-18
            f64::from_bits(0xBBD89183964103F4), // -2.081025621058701e-20
        ],
    },
    // [24.347343065320896, 27.48893571891069)
    ZeroPoly {
        c: f64::from_bits(0x4039EC46F3E80146), // 2.5922957653180923e1
        k_hi: [
            F64x2::from_bits(0x3C845D52C7C075B2, 0xB904F195416FF08D), // 3.532685320807229012585113578800e-17
            F64x2::from_bits(0x3FC40F8FFDF09A5F, 0x3C58EB882BA0A2C7), // 1.567249288502407854242155404552e-1
            F64x2::from_bits(0xBF68C37A29C45897, 0xBC026DE49F51CE5E), // -3.022898292452568354946968872704e-3
            F64x2::from_bits(0xBF9AAB099314B209, 0xBC28357573FBE6BD), // -2.604308090410437505739041758308e-2
            F64x2::from_bits(0x3F406F735CC0F6B4, 0xBBDC17222E5E7CE9), // 5.015672009062916510899211023265e-4
            F64x2::from_bits(0x3F552D29A06802E0, 0xBBE53A6F72A8424E), // 1.292505885089827130771934694410e-3
            F64x2::from_bits(0xBEF9FC04C675C0F9, 0x3B9904FA38A8C80F), // -2.478070055417309480136378469685e-5
            F64x2::from_bits(0xBEFFE48825ED3C44, 0x3B8DAE6298AC5539), // -3.041525121299746482681102116932e-5
            F64x2::from_bits(0x3EA36BD2D58FF354, 0xBAEC02B78CC58AE9), // 5.787964411082431619432863160275e-7
            F64x2::from_bits(0x3E9BE87E2CACFD52, 0x3B3022969091D843), // 4.158642183289454907685788931443e-7
            F64x2::from_bits(0xBE40D1FDEDC6461F, 0x3AB0306443D8D6D8), // -7.832553654954706487226531968176e-9
        ],
        k_lo: [
            f64::from_bits(0xBE2FDBEB105A9406), // -3.708882203813804e-9
            f64::from_bits(0x3DD2F44D00E4B62C), // 6.89553416305312e-11
            f64::from_bits(0x3DB9921363EFEF1B), // 2.3256332837005183e-11
            f64::from_bits(0xBD5DF6B803480AF7), // -4.258104313631636e-13
            f64::from_bits(0xBD3E6BF760BF5F76), // -1.0807974404337205e-13
            f64::from_bits(0x3CE184804994E28B), // 1.944842343926561e-15
            f64::from_bits(0x3CBBE2443CC7FEE7), // 3.869662051972376e-16
            f64::from_bits(0xBC5F7ECB00E83D35), // -6.82945130314664e-18
            f64::from_bits(0xBC33BDC188864A9A), // -1.070174457611914e-18
            f64::from_bits(0x3BD596C48BE6CA5E), // 1.828654056372008e-20
        ],
    },
    // [27.48893571891069, 30.630528372500482)
    ZeroPoly {
        c: f64::from_bits(0x403D106449616C4F), // 2.9064030252728397e1
        k_hi: [
            F64x2::from_bits(0x3CA3BC249963D4C1, 0xB941925C39C8D9F0), // 1.369386094737845000728640104157e-16
            F64x2::from_bits(0xBFC2F206E49909C7, 0xBC66277503C2B1D5), // -1.480110756111356714081547672415e-1
            F64x2::from_bits(0x3F64DBF6A9FB803A, 0xBBB65D633153EFF5), // 2.546293035138117057929745178975e-3
            F64x2::from_bits(0x3F99336443318ED1, 0x3C23E7228A486375), // 2.461010608784769805413842494336e-2
            F64x2::from_bits(0xBF3BB6AA3D4E9E10, 0x3BC19C728A4FBC8A), // -4.228749870045242770587477183945e-4
            F64x2::from_bits(0xBF540AAA5D94BD8D, 0xBBD852B2412167A5), // -1.223245976074458016413238984454e-3
            F64x2::from_bits(0x3EF5F61B6661296D, 0xBB8AE2F748D780AE), // 2.094398077926702767045986093148e-5
            F64x2::from_bits(0x3EFE4158391F2BE4, 0x3B79CFC6747FBE6D), // 2.885365717712394790614967165158e-5
            F64x2::from_bits(0xBEA07A7A0745F63C, 0xBB439FA9C3F5B795), // -4.910953195627749797296985147844e-7
            F64x2::from_bits(0xBE9A8EA97B66D1F9, 0x3B3CCC1B1B7F7532), // -3.957342146435761890641971899544e-7
            F64x2::from_bits(0x3E3CB10B3AFF012F, 0x3AB51FA3BE9E2F63), // 6.680278483108145969404469598637e-9
        ],
        k_lo: [
            f64::from_bits(0x3E2E6D78877D34BC), // 3.5422412190197893e-9
            f64::from_bits(0xBDD0444478C004D0), // -5.917779997588326e-11
            f64::from_bits(0xBDB8850C8A24FDC2), // -2.2300557784112628e-11
            f64::from_bits(0x3D59E4FA5F26E2FB), // 3.6798220106619625e-13
            f64::from_bits(0x3D3D4B235AAA17C8), // 1.0407144732354132e-13
            f64::from_bits(0xBCDE806889119426), // -1.6931786575020533e-15
            f64::from_bits(0xBCBAF647B531CF1F), // -3.74173354440194e-16
            f64::from_bits(0x3C5BA099236D7E94), // 5.990723695409749e-18
            f64::from_bits(0x3C3329FAEFD9FCD3), // 1.0388817217257643e-18
            f64::from_bits(0xBBD314263C6A7A61), // -1.616029463554962e-20
        ],
    },
    // [30.630528372500482, 33.772121026090275)
    ZeroPoly {
        c: f64::from_bits(0x40401A4420E4ABEE), // 3.220520411649328e1
        k_hi: [
            F64x2::from_bits(0xBCC069EACEA8FF0B, 0xB96FD1B1319B777A), // -4.555727774065143303849023423341e-16
            F64x2::from_bits(0x3FC1FF5EBDDD3C3A, 0x3C42F0FF4F7EA2BA), // 1.406057765075006703023260142885e-1
            F64x2::from_bits(0xBF61E20353246230, 0xBBE8D9079CC73F31), // -2.182966703127918504006790036943e-3
            F64x2::from_bits(0xBF97F3506D4A1233, 0x3C2E26B27420FB96), // -2.338910738991923769934686943088e-2
            F64x2::from_bits(0x3F37C65C9302C3DF, 0xBBCA88987A782B16), // 3.627754222277883339732515733510e-4
            F64x2::from_bits(0x3F53117816335152, 0xBBD871B513775DDC), // 1.163832924068481929875206048431e-3
            F64x2::from_bits(0xBEF2DF9AFA5211DA, 0x3B92652CE0F0818C), // -1.799913265589253092541318712305e-5
            F64x2::from_bits(0xBEFCD5D4A9D789D1, 0xBB9648D598A0F3FF), // -2.749946235337271686071216305107e-5
            F64x2::from_bits(0x3E9C672D4D6F74F0, 0xBB3E7BB566D2D658), // 4.232382030936481653926769203750e-7
            F64x2::from_bits(0x3E995E4B102167CE, 0x3B25E966A23F77BE), // 3.780176173244688744589776979342e-7
            F64x2::from_bits(0xBE38D22636133E98, 0xBAD462140E947772), // -5.779065088628769134170549987333e-9
        ],
        k_lo: [
            f64::from_bits(0xBE2D24FC943744DA), // -3.392863907990917e-9
            f64::from_bits(0x3DCC462D9D2DBC67), // 5.143034925656279e-11
            f64::from_bits(0x3DB78EAFA2309A56), // 2.1425300888363044e-11
            f64::from_bits(0xBD569F5D1E47D78C), // -3.214852685893159e-13
            f64::from_bits(0xBD3C3C45D621A935), // -1.0031243611648434e-13
            f64::from_bits(0x3CDACC9BDFD31266), // 1.4876574113300829e-15
            f64::from_bits(0x3CBA1335B01269F3), // 3.6186384389687444e-16
            f64::from_bits(0xBC586B11C5C39F80), // -5.294861757610048e-18
            f64::from_bits(0xBC3297F5A9EB9FAC), // -1.00796065728832e-18
            f64::from_bits(0x3BD0F69E0154037A), // 1.436851504801408e-20
        ],
    },
    // [33.772121026090275, 36.91371367968007)
    ZeroPoly {
        c: f64::from_bits(0x4041AC588C944279), // 3.534645230521432e1
        k_hi: [
            F64x2::from_bits(0xBC7BC4FD81E2B2C5, 0x391BF56A50DB08B2), // -2.408619589837893008575852386589e-17
            F64x2::from_bits(0xBFC12DD55D4BE2B3, 0xBC48318D12A5D220), // -1.342112260388340387458641428592e-1
            F64x2::from_bits(0x3F5F1AEE31818D50, 0x3BD7C71D2496CB13), // 1.898510561681405354340649200507e-3
            F64x2::from_bits(0x3F96DE64242A8310, 0xBC35DD175356FA59), // 2.233273001423002771796920843267e-2
            F64x2::from_bits(0xBF34AFDF89FCA62D, 0x3BD0B74F3D7712C9), // -3.156586408175923003230786957920e-4
            F64x2::from_bits(0xBF5238CFC13AC771, 0x3BF8D99C3C08C4A5), // -1.112177739667843384103957583008e-3
            F64x2::from_bits(0x3EF0719D13E00E5A, 0x3B9E352E5043D18F), // 1.568203264669244933873261830463e-5
            F64x2::from_bits(0x3EFB974781A52670, 0xBB99B5246EE362B1), // 2.631276593348944378263545149688e-5
            F64x2::from_bits(0xBE98CC82A70D73C2, 0x3B1C586D4B7A4C1C), // -3.695319384289664265740001493850e-7
            F64x2::from_bits(0xBE9850AE878BF9CD, 0x3B1222FA14C1757A), // -3.623241649346734789158640159681e-7
            F64x2::from_bits(0x3E35BBA73E277650, 0xBA6446466A279A10), // 5.060113191709251436365420730276e-9
        ],
        k_lo: [
            f64::from_bits(0x3E2BFE1396C7A227), // 3.2587543134775104e-9
            f64::from_bits(0xBDC8D7DBE436BC48), // -4.51895157810355e-11
            f64::from_bits(0xBDB6ADFBC8B313E6), // -2.0626997110119578e-11
            f64::from_bits(0x3D53F51DA4991F27), // 2.8361289943313717e-13
            f64::from_bits(0x3D3B40CAC1AC7E77), // 9.682243920127915e-14
            f64::from_bits(0xBCD7BF432242C38B), // -1.318229865874808e-15
            f64::from_bits(0xBCB93C43853ADA5B), // -3.5021159971987005e-16
            f64::from_bits(0x3C55BD42CDA7F9C5), // 4.713959385195399e-18
            f64::from_bits(0x3C320AEBE6B32667), // 9.780946704803123e-19
            f64::from_bits(0xBBCE5997A0637F2B), // -1.2853712221256787e-20
        ],
    },
    // [36.91371367968007, 40.05530633326986)
    ZeroPoly {
        c: f64::from_bits(0x40433E6ECF5CB221), // 3.848775665308154e1
        k_hi: [
            F64x2::from_bits(0x3CBA82B3A8997D82, 0x393AB2F3287D39B6), // 3.679078413110179927260385612687e-16
            F64x2::from_bits(0x3FC0768257DAD569, 0x3C65C008F799FF5F), // 1.286166123369796671040390363835e-1
            F64x2::from_bits(0xBF5B602A7BEAAD97, 0x3BDD5A9E5AE0B919), // -1.670876968698253272940447659196e-3
            F64x2::from_bits(0xBF95EBC22EFD092B, 0x3BD8EE8893D9F668), // -2.140715991777956612282914670186e-2
            F64x2::from_bits(0x3F3236A604142F7A, 0x3BDB7383389AF036), // 2.779155075394096096884294995823e-4
            F64x2::from_bits(0x3F517A482FAA8D84, 0xBBD102BA8F113DEC), // 1.066751951534467855200346565982e-3
            F64x2::from_bits(0xBEECFC00890788D6, 0x3B8A8E716399FDCF), // -1.382083090190439167980374143684e-5
            F64x2::from_bits(0xBEFA7D7B128AC4F3, 0xBB6871EF82550861), // -2.526298445343784973321760663021e-5
            F64x2::from_bits(0x3E95E419F1B56F9B, 0xBB365F5AE446791B), // 3.262016307814389945728462981339e-7
            F64x2::from_bits(0x3E97603CEA232781, 0xBB27AD9BEBBB4270), // 3.483284932824317798841658934313e-7
            F64x2::from_bits(0xBE333A0B436BE2D3, 0x3ACC98AA169D0B7A), // -4.476572937637116050723662718596e-9
        ],
        k_lo: [
            f64::from_bits(0xBE2AF4641764AACB), // -3.137934519123691e-9
            f64::from_bits(0x3DC609AC8C77806A), // 4.008650491267446e-11
            f64::from_bits(0x3DB5E11B79522312), // 1.9899130592048745e-11
            f64::from_bits(0xBD51C27542D9A321), // -2.523791204524099e-13
            f64::from_bits(0xBD3A585C7291378F), // -9.359681256645383e-14
            f64::from_bits(0x3CD534C2C0BACA7E), // 1.1771748405309583e-15
            f64::from_bits(0x3CB8729B831647DF), // 3.3927979073007366e-16
            f64::from_bits(0xBC537C8A9ED22FD7), // -4.225458997316442e-18
            f64::from_bits(0xBC3184841F4F6E3A), // -9.496332228818848e-19
            f64::from_bits(0x3BCB4FA5E446F94F), // 1.156671137475253e-20
        ],
    },
    // [40.05530633326986, 43.19689898685966)
    ZeroPoly {
        c: f64::from_bits(0x4044D0867EC213F3), // 4.162910446621381e1
        k_hi: [
            F64x2::from_bits(0xBCB71414F52D813B, 0x395EC360D6E463F7), // -3.202777597232382929834738547878e-16
            F64x2::from_bits(0xBFBFA8B3F9AE4374, 0xBC5B80DAEB480DF4), // -1.236679539272367904476037421988e-1
            F64x2::from_bits(0x3F5856073B7FA5AF, 0xBBF9CC6931E23BCF), // 1.485354483515482050211871837665e-3
            F64x2::from_bits(0x3F9514E652EB2E96, 0xBC38E4CB2C3A351C), // 2.058753853888946886928187981165e-2
            F64x2::from_bits(0xBF3032298718EB0F, 0xBBD103147CCFE20B), // -2.471305261654713152375261245636e-4
            F64x2::from_bits(0xBF50D1153FDE4319, 0xBBFCC81D970A2C4B), // -1.026411773157492846761143426534e-3
            F64x2::from_bits(0x3EE9CB455C29D7A9, 0xBB69E569AFE43CAD), // 1.229955058089009744446652482114e-5
            F64x2::from_bits(0x3EF982167B1C5280, 0x3B59A49909745452), // 2.432647278912064914438093295014e-5
            F64x2::from_bits(0xBE938191EF68EA15, 0xBB28213208E7711F), // -2.906640327625192373290182489599e-7
            F64x2::from_bits(0xBE968865345B078E, 0xBB2E744E19DEF81B), // -3.357647993911996660582203738964e-7
            F64x2::from_bits(0x3E312960010D0E20, 0xBAD0865248B25093), // 3.995751299987148688698359213315e-9
        ],
        k_lo: [
            f64::from_bits(0x3E2A040367808A8B), // 3.028623404164557e-9
            f64::from_bits(0xBDC3B67CC0621F20), // -3.585744899606881e-11
            f64::from_bits(0xBDB526037D9B7AF1), // -1.9234440299577135e-11
            f64::from_bits(0x3D4FDA06A97C6611), // 2.2631968587122608e-13
            f64::from_bits(0x3D3981D68958EE6E), // 9.061970664490693e-14
            f64::from_bits(0xBCD311EE641CDFAE), // -1.0586000858672254e-15
            f64::from_bits(0xBCB7B64970D13DD2), // -3.290709110209623e-16
            f64::from_bits(0x3C5193CC012C0DF4), // 3.811476224225822e-18
            f64::from_bits(0x3C31056F6A38F3F6), // 9.2272279806806e-19
            f64::from_bits(0xBBC8B6D1F843BA07), // -1.046684600794335e-20
        ],
    },
    // [43.19689898685966, 46.33849164044945)
    ZeroPoly {
        c: f64::from_bits(0x4046629F4E1E0321), // 4.477048660722199e1
        k_hi: [
            F64x2::from_bits(0xBCBAE2D399EAFCC6, 0x3946B712A4F82541), // -3.731187758581486799338405727576e-16
            F64x2::from_bits(0x3FBE8727C572A2C3, 0xBC27499384D0D237), // 1.192498070764811812021874814500e-1
            F64x2::from_bits(0xBF55D1EF092AB039, 0xBBC325012256146E), // -1.331790383726007154358634501742e-3
            F64x2::from_bits(0xBF9454E7A7395637, 0x3C351B3A50E28076), // -1.985513647226452387741356111292e-2
            F64x2::from_bits(0x3F2D0CC3A7FA6B02, 0xBBC8D5B2311F4319), // 2.216328466502979092682896931067e-4
            F64x2::from_bits(0x3F50398D2CBD02E0, 0xBBF29043E49A001E), // 9.902838382488884298433352487273e-4
            F64x2::from_bits(0xBEE725B0909F9B34, 0x3B8CCB9C6BFD8119), // -1.103745718817952862674619453971e-5
            F64x2::from_bits(0xBEF8A022F5FEE40A, 0xBB949F8C14A7127A), // -2.348473878985641472124703697401e-5
            F64x2::from_bits(0x3E9185B1EA97A389, 0xBB1E099ABDEFCC9F), // 2.611018126973077814943977273071e-7
            F64x2::from_bits(0x3E95C5711673D924, 0x3B3D17B79029F450), // 3.244170074379334068556445732328e-7
            F64x2::from_bits(0xBE2EDFF6514D0C70, 0x3ACED865FA685394), // -3.594305861977030464457748175477e-9
        ],
        k_lo: [
            f64::from_bits(0xBE29298DA3A6E8BA), // -2.929279288937857e-9
            f64::from_bits(0x3DC1C347059B8741), // 3.2310349448850364e-11
            f64::from_bits(0x3DB47AB8DEB5230D), // 1.8625890690230213e-11
            f64::from_bits(0xBD4CC249F5F1276F), // -2.0434456651861113e-13
            f64::from_bits(0xBD38BBC49E2CCE02), // -8.78709332874267e-14
            f64::from_bits(0x3CD142014A973DAF), // 9.580021334413438e-16
            f64::from_bits(0x3CB706C48740A64A), // 3.1955599672417033e-16
            f64::from_bits(0xBC4FE4D0B874C417), // -3.4579337898465942e-18
            f64::from_bits(0xBC308DCEF05A25D1), // -8.97390825204952e-19
            f64::from_bits(0x3BC67B1F8439CE33), // 9.521052521138474e-21
        ],
    },
    // [46.33849164044945, 49.480084294039244)
    ZeroPoly {
        c: f64::from_bits(0x4047F4B904DC9A53), // 4.791189633151648e1
        k_hi: [
            F64x2::from_bits(0xBC78E44442565953, 0xB91E0B0008451F03), // -2.159007916168211903744625549971e-17
            F64x2::from_bits(0xBFBD82939AB62339, 0x3C577133F0CA2735), // -1.152736904824792717699172803682e-1
            F64x2::from_bits(0x3F53B5A548456741, 0xBBFB9ED457AABF9A), // 1.202975662712949102030762029453e-3
            F64x2::from_bits(0x3F93A7FF1622DEF8, 0xBC35F259372B469E), // 1.919554302796108605598253955057e-2
            F64x2::from_bits(0xBF2A3EBC476A608E, 0xBBBB0649B931E292), // -2.002339203570662994906572910559e-4
            F64x2::from_bits(0xBF4F61ADDE3A8C61, 0x3BEC054F8EEB22E4), // -9.576891865957886439408696579062e-4
            F64x2::from_bits(0x3EE4EC45DA0478BF, 0x3B8321D65040F472), // 9.976835656785495330975660858243e-6
            F64x2::from_bits(0x3EF7D3B28159F1FE, 0x3B9ED54382A127BB), // 2.272314312400276244598908034812e-5
            F64x2::from_bits(0xBE8FB3F21BCD7DDD, 0x3AD2BB3342A20A73), // -2.362051087180024969543608960592e-7
            F64x2::from_bits(0xBE95145A876CFAEC, 0x3B3CB7E4483DCCD3), // -3.141091222020475523775268699623e-7
            F64x2::from_bits(0x3E2BF64CDF7B0CB9, 0x3AC5E497EC50BE5C), // 3.255218091473714039345466460319e-9
        ],
        k_lo: [
            f64::from_bits(0x3E2862204AA6D7B4), // 2.838590325652596e-9
            f64::from_bits(0xBDC01C0D48400E30), // -2.93031510786137e-11
            f64::from_bits(0xBDB3DD6CA4860028), // -1.806705677661734e-11
            f64::from_bits(0x3D4A2020949219C1), // 1.856328220864194e-13
            f64::from_bits(0x3D3804A2619D34C9), // 8.53294421506273e-14
            f64::from_bits(0xBCCF6A07C9F2AE22), // -8.719186858230484e-16
            f64::from_bits(0xBCB6633DF2004F08), // -3.106912499568727e-16
            f64::from_bits(0x3C4D169863B792B4), // 3.153755769940517e-18
            f64::from_bits(0x3C301D72CAE8375C), // 8.735976810740115e-19
            f64::from_bits(0xBBC48C6D482752F8), // -8.70264626146141e-21
        ],
    },
];

// GENERATE: bessel_table y1 15 2 1 11 10
// Y1_TBL[i] has the coefficients `K0, ..., K20` such that
// f(c + d) ~= K0 + K1 * d + ... + K20 * d^20, where `c` is very close to the
// (i + 1)-th positive zero of y1, for the values of the interval of that
// zero, [(i + 0.25) * π, (i + 1.25) * π) (the first one starts at 1)
// f(x) = x * y1(x) - 2/π * x * j1(x) * ln(x / c) in the first 2 intervals
// f(x) = y1(x) in the other ones
static Y1_TBL: [ZeroPoly; 15] = [
    // [1, 3.9269908169872414)
    ZeroPoly {
        c: f64::from_bits(0x400193BED4DFF243), // 2.197141326031017e0
        k_hi: [
            F64x2::from_bits(0x3C8FD52A9B5D6777, 0x392EC94FF9CE11B0), // 5.522089969249562975687058459197e-17
            F64x2::from_bits(0x3FE9480EF7F8CC11, 0x3C6F6306FB94E433), // 7.900462001038041218286862933344e-1
            F64x2::from_bits(0x3FD1433054BB1EF2, 0xBC67167A7839E9D1), // 2.697258784723849289536597863674e-1
            F64x2::from_bits(0xBFA37D71C727F092, 0x3C41B14187782B16), // -3.806644016364758661666418939056e-2
            F64x2::from_bits(0xBFA35F4193C6875E, 0x3C42A4F9C71B3C10), // -3.783612183475669263816228646649e-2
            F64x2::from_bits(0xBF65DFD7833147B2, 0xBC06B3AFCB156E1E), // -2.670212672037452683383893894839e-3
            F64x2::from_bits(0x3F5C75A58B296355, 0xBBF2F7EE1589C186), // 1.737033523246907566667588692309e-3
            F64x2::from_bits(0x3F25D972CD4641B2, 0x3BC03B8E9B6CEC87), // 1.666977538090759235427235727068e-4
            F64x2::from_bits(0xBF04847D3C4C4B48, 0x3BA3168D3B95BFE9), // -3.913409412894267040166367207011e-5
            F64x2::from_bits(0xBECED16A66CFD494, 0x3B3A08277FC78994), // -3.673795438969573238368370816008e-6
            F64x2::from_bits(0x3EA180866E42112B, 0x3B1D4BFA269DDF38), // 5.216017737873714314111996217730e-7
        ],
        k_lo: [
            f64::from_bits(0x3E6840F1D7D33D9F), // 4.517601846222952e-8
            f64::from_bits(0xBE33B348CDF48D0B), // -4.586840434629385e-9
            f64::from_bits(0xBDF8DFB3B9133BC4), // -3.6196195439813323e-10
            f64::from_bits(0x3DBF6E53DE065F02), // 2.8586298152010132e-11
            f64::from_bits(0x3D82147612456666), // 2.055449684366002e-12
            f64::from_bits(0xBD42B13B6AB213D0), // -1.3281687130095255e-13
            f64::from_bits(0xBD03ACD300535A1B), // -8.737497225252558e-15
            f64::from_bits(0x3CC146E3717276CF), // 4.795305265951287e-16
            f64::from_bits(0x3C80A35D4C07AEEB), // 2.886257612967055e-17
            f64::from_bits(0xBC3ECCAEEDE2AF60), // -1.6696466371142793e-18
        ],
    },
    // [3.9269908169872414, 7.0685834705770345)
    ZeroPoly {
        c: f64::from_bits(0x4015B7FE4E87B02E), // 5.429681040794135e0
        k_hi: [
            F64x2::from_bits(0x3CCBB638A8D0BB1D, 0xB968DCE933F8CFA8), // 7.691570173653176440805592736950e-16
            F64x2::from_bits(0xBFFA0E485A4AC487, 0xBC5424D57228D6D2), // -1.628486969664921128957804018235e0
            F64x2::from_bits(0xBFC5D783DB512373, 0xBC47E37700689729), // -1.706394978412180988645737757259e-1
            F64x2::from_bits(0x3FC9B3495289B112, 0xBC21A373CDB1F1FF), // 2.007838871971974743162158216529e-1
            F64x2::from_bits(0x3F9CD42ACFDB1CA7, 0xBC27FDE776EFB137), // 2.815310376788873335447679903558e-2
            F64x2::from_bits(0xBF7CE25038511BC4, 0x3BFB735014B174E0), // -7.051766738249424283622061236152e-3
            F64x2::from_bits(0xBF55F3E64EF796FC, 0xBB71DAC19F11FA3F), // -1.339888487444650632077583405157e-3
            F64x2::from_bits(0x3F1A014A186C52DA, 0x3BA311ED45015D4A), // 9.920134298753037532966555214215e-5
            F64x2::from_bits(0x3F000C0679F13899, 0xBB9489545972D326), // 3.060717357837687042491825215296e-5
            F64x2::from_bits(0xBE9DF4A337503F54, 0x3B2EDD660167F9FB), // -4.463734549217924084705609313572e-7
            F64x2::from_bits(0xBE9B9AFE813F947F, 0xBB08C5B7353600B0), // -4.113531997237270795997464648394e-7
        ],
        k_lo: [
            f64::from_bits(0xBE3210520609A227), // -4.205794907308637e-9
            f64::from_bits(0x3E2F46E8B88DCBC0), // 3.6411206865836364e-9
            f64::from_bits(0x3DD59DCFDC4F4DD5), // 7.864019771958627e-11
            f64::from_bits(0xBDB9188C6CB2C505), // -2.282458145545518e-11
            f64::from_bits(0xBD65A036CE34B82A), // -6.146432345258543e-13
            f64::from_bits(0x3D3E02794817BED8), // 1.0661574062900905e-13
            f64::from_bits(0x3CEC1AB9D7C489F3), // 3.120215001565836e-15
            f64::from_bits(0xBCBBC8B4A3E4522D), // -3.855805445866566e-16
            f64::from_bits(0xBC69D8913D2F4536), // -1.1208901069165592e-17
            f64::from_bits(0x3C348AB90E4BED1D), // 1.1135778837738635e-18
        ],
    },
    // [7.0685834705770345, 10.210176124166829)
    ZeroPoly {
        c: f64::from_bits(0x40213127AE6169B4), // 8.596005868331169e0
        k_hi: [
            F64x2::from_bits(0xBC963BC010B45F46, 0x393EA23630C8635A), // -7.713759989498053960149631972772e-17
            F64x2::from_bits(0x3FD15F993FCEAB5C, 0x3C7346BC5A82E478), // 2.714598773115335540002872563907e-1
            F64x2::from_bits(0xBF902B3933CF21A6, 0xBC3A82CA93B4694C), // -1.578988436429686897319614605980e-2
            F64x2::from_bits(0xBFA6395DFE49FCD4, 0xBC2433C9A5FEB405), // -4.340642670740071983683806433950e-2
            F64x2::from_bits(0x3F63CED2A2E6917A, 0xBC0F678499AA98BE), // 2.417956732829459847803293846847e-3
            F64x2::from_bits(0x3F607A678D6000BC, 0xBC019EEF1FF7CE01), // 2.011492014389741156355995805622e-3
            F64x2::from_bits(0xBF1B50D7E1D32021, 0xBB949FD8C90926DB), // -1.042014850610732999041989147909e-4
            F64x2::from_bits(0xBF06F7BAB104F31D, 0xBB91D2A6E3B88761), // -4.380739676654850684102501267807e-5
            F64x2::from_bits(0x3EC176E72BF94ECD, 0x3B613E3BB4964491), // 2.081926453270871992330206087733e-6
            F64x2::from_bits(0x3EA2BECB2B6B9786, 0x3B4E1E0C45D47A39), // 5.586531049677760347191018503931e-7
            F64x2::from_bits(0xBE5A384EEC01045A, 0xBAFADC97A2537C63), // -2.441923530071958140470350166420e-8
        ],
        k_lo: [
            f64::from_bits(0xBE341E7A9220550C), // -4.684333173109147e-9
            f64::from_bits(0x3DE9E32850E84A75), // 1.8835514437537138e-10
            f64::from_bits(0x3DBEC40B4870B79E), // 2.798132951689526e-11
            f64::from_bits(0xBD72686932668841), // -1.0463654221528287e-12
            f64::from_bits(0xBD416F737DB1D4D2), // -1.2388565553129188e-13
            f64::from_bits(0x3CF3189065E76DFE), // 4.240153414953438e-15
            f64::from_bits(0x3CBFA7820C2B18C5), // 4.3929204881870837e-16
            f64::from_bits(0xBC7143DA9EB0534A), // -1.4975047785660767e-17
            f64::from_bits(0xBC32CD1358BDAA88), // -1.0192083969526342e-18
            f64::from_bits(0x3BDA86A114D03CFA), // 2.2468307363978386e-20
        ],
    },
    // [10.210176124166829, 13.351768777756622)
    ZeroPoly {
        c: f64::from_bits(0x40277F9138D43206), // 1.1749154830839881e1
        k_hi: [
            F64x2::from_bits(0x3C5F96D32C02F147, 0xB8EF1284C1B1C7AD), // 6.849807120770496002079016654612e-18
            F64x2::from_bits(0xBFCDC14EA14E89F9, 0x3C5221104CBAFD09), // -2.324617660170387417152831022975e-1
            F64x2::from_bits(0x3F84429FEF5B5FBB, 0x3C27A6F8F4FF8B98), // 9.892701618284030790015009938592e-3
            F64x2::from_bits(0x3FA367D7D608E4BB, 0xBC4F3801F21E0AA5), // 3.790163505295525781112591744195e-2
            F64x2::from_bits(0xBF59D6EB2BC49E34, 0xBBDFBE2B7AF5F77E), // -1.577119502209967306058794711249e-3
            F64x2::from_bits(0xBF5DC4F991B3DB86, 0x3BE326DA984FB445), // -1.816982002137841492763429907460e-3
            F64x2::from_bits(0x3F1315EC04D6E6C0, 0xBBB71F69A492ECA6), // 7.280591054028683266585040664520e-5
            F64x2::from_bits(0x3F0571814A1AA2DE, 0xBBA57FFB40A620FC), // 4.089999971529864933880813260003e-5
            F64x2::from_bits(0xBEBA2977FA42F557, 0x3B510D11D0A3B3C1), // -1.559375939527129421848802691678e-6
            F64x2::from_bits(0xBEA1E86423084C15, 0x3B1D14DEE75CEEF9), // -5.336933722198205676803351710591e-7
            F64x2::from_bits(0x3E54A7B82D436A7C, 0x3AD41DB4568D66D6), // 1.923661126298179409600564201875e-8
        ],
        k_lo: [
            f64::from_bits(0x3E3390660A0A68BD), // 4.555111982616221e-9
            f64::from_bits(0xBDE549E8B4B2DBBB), // -1.5489573022599694e-10
            f64::from_bits(0xBDBE32CF67D2B021), // -2.7465355067460198e-11
            f64::from_bits(0x3D6EFF5953223D82), // 8.810007084714545e-13
            f64::from_bits(0x3D415E5784349ED2), // 1.2341077773721443e-13
            f64::from_bits(0xBCF0D696D6AA5622), // -3.738840151570978e-15
            f64::from_bits(0xBCBEF0D0CA9515BB), // -4.2938827358971996e-16
            f64::from_bits(0x3C6C69CF5C6FA2D7), // 1.2322313451272066e-17
            f64::from_bits(0x3C356AAE4C695CE1), // 1.1610028305377128e-18
            f64::from_bits(0xBBE3AB71906AF004), // -3.332177063763767e-20
        ],
    },
    // [13.351768777756622, 16.493361431346415)
    ZeroPoly {
        c: f64::from_bits(0x402DCB7D88DE848B), // 1.4897442128336726e1
        k_hi: [
            F64x2::from_bits(0x3CA213201464272B, 0x3920733C1B529666), // 1.254184828676617677430691554314e-16
            F64x2::from_bits(0x3FCA7022BE084D99, 0xBC4E3247E4F1BDE5), // 2.065471103565925914970290242458e-1
            F64x2::from_bits(0xBF7C650B6B8310E2, 0x3C1434BF144E3F2E), // -6.932301148655480100662352161631e-3
            F64x2::from_bits(0xBFA163191C30AA62, 0xBC1EBD3877DB2544), // -3.395918339841298747043127660385e-2
            F64x2::from_bits(0x3F526B045287DDE1, 0x3BEA3C770EEF2938), // 1.124147625999210015452502402579e-3
            F64x2::from_bits(0x3F5B17602840ABF4, 0xBBFFBD063AAA2F11), // 1.653522399474691217271707013477e-3
            F64x2::from_bits(0xBF0C0A9CEE3C8448, 0x3B793B59E81AA624), // -5.348483480661269013049843379632e-5
            F64x2::from_bits(0xBF03E398CBC472C4, 0xBBAFDE41AC782C6E), // -3.793535276683726512408899956764e-5
            F64x2::from_bits(0x3EB3F35DB1FF1E1D, 0x3B59104D61B11FAA), // 1.189151312518366643008054848817e-6
            F64x2::from_bits(0x3EA0E9B612DBBBB6, 0xBB4BE037C28ADEFC), // 5.040447257057908914366128126292e-7
            F64x2::from_bits(0xBE5056BABCDA86C5, 0xBAD2555D93E02B75), // -1.521668107371637180358286089614e-8
        ],
        k_lo: [
            f64::from_bits(0xBE32C1A8C8E79CD3), // -4.3670837072231724e-9
            f64::from_bits(0x3DE161B6AB3A3F1C), // 1.2646846557301765e-10
            f64::from_bits(0x3DBD4CAA0DC7614A), // 2.6647712561882862e-11
            f64::from_bits(0xBD69FDDA74BDA234), // -7.387261184790911e-13
            f64::from_bits(0xBD4101BE84786246), // -1.2084067659944283e-13
            f64::from_bits(0x3CECDC495EF2DBD1), // 3.2041585552570316e-15
            f64::from_bits(0x3CBE83D8616D2C3B), // 4.2348097883264117e-16
            f64::from_bits(0xBC68CCAD8AB5FD6A), // -1.0755078701901046e-17
            f64::from_bits(0xBC3539B1A166E459), // -1.1506294330943118e-18
            f64::from_bits(0x3BE14E72DD72720F), // 2.9318250100249615e-20
        ],
    },
    // [16.493361431346415, 19.634954084936208)
    ZeroPoly {
        c: f64::from_bits(0x40320B1C695F1E3B), // 1.8043402276727857e1
        k_hi: [
            F64x2::from_bits(0xBCB39D4C41D5839F, 0xB958FF9BB96B5AA0), // -2.722051034851249658697830812014e-16
            F64x2::from_bits(0xBFC80781C32422E7, 0x3C5DDBCC2B07582A), // -1.877290919149096839116333995324e-1
            F64x2::from_bits(0x3F754EDA697A0133, 0x3BE0CA5E235714DD), // 5.202153369851014585301416544774e-3
            F64x2::from_bits(0x3F9FBE6DF840847E, 0xBC318DC9061C900F), // 3.099986865720837691293541562612e-2
            F64x2::from_bits(0xBF4BE318D6127747, 0x3BE670D8F040C5C9), // -8.510466833688633540597646123560e-4
            F64x2::from_bits(0xBF58EFEE4094379A, 0xBBF2D13995A104D5), // -1.522047680114968643606317159402e-3
            F64x2::from_bits(0x3F059145B4F0E518, 0xBBA1014BA37852B8), // 4.113668421052758246441845512806e-5
            F64x2::from_bits(0x3F0282D26A74C367, 0x3BA43F4CF104E5AF), // 3.530697477527912974464208293721e-5
            F64x2::from_bits(0xBEAF56C29D9ADDA8, 0xBB3EC7546700925B), // -9.339722126443645327774447490070e-7
            F64x2::from_bits(0xBE9FDD03174F603E, 0xBB3F68725ED74E19), // -4.748005928956742835902025005345e-7
            F64x2::from_bits(0x3E4A44A79081C812, 0xBAC9DB2CF520DABB), // 1.223207536668488789255440098867e-8
        ],
        k_lo: [
            f64::from_bits(0x3E31DF653318DAEE), // 4.1612977938072865e-9
            f64::from_bits(0xBDDC96E9D0301C9A), // -1.0400801404764774e-10
            f64::from_bits(0xBDBC3439E3BAEEB4), // -2.5651396140123802e-11
            f64::from_bits(0x3D65D38B33B1FE59), // 6.203419730740883e-13
            f64::from_bits(0x3D40833A0A0EE7AE), // 1.1732911075968682e-13
            f64::from_bits(0xBCE8AFD3A3A999A8), // -2.740787942105684e-15
            f64::from_bits(0xBCBDD71AE45C5249), // -4.1411671439048185e-16
            f64::from_bits(0x3C658DD5956EAD95), // 9.347574918503095e-18
            f64::from_bits(0x3C34DF75F9F38856), // 1.1315218470628061e-18
            f64::from_bits(0xBBDE6325775DC9FC), // -2.5739036180452915e-20
        ],
    },
    // [19.634954084936208, 22.776546738526)
    ZeroPoly {
        c: f64::from_bits(0x40353025492188CD), // 2.1188068934142212e1
        k_hi: [
            F64x2::from_bits(0xBCAB20138DD047EC, 0xB8F3AA7E1213624B), // -1.882195674994624986251799025106e-16
            F64x2::from_bits(0x3FC62D94D97E859C, 0xBC516689EDCBC66B), // 1.732660352691198741864127171373e-1
            F64x2::from_bits(0xBF70BF61480702D0, 0xBC0A8B429BE98B11), // -4.088764195728969610423439938818e-3
            F64x2::from_bits(0xBF9D5F857A2A6108, 0xBC29A7573DBE3AFC), // -2.868469770564366430755145860054e-2
            F64x2::from_bits(0x3F46081B0B7FE52B, 0x3BECCDD29882E929), // 6.723529868387970100608951398608e-4
            F64x2::from_bits(0x3F57307B03E248F9, 0x3BEAD36BF9300007), // 1.415367252353193132369839089678e-3
            F64x2::from_bits(0xBF0132C0AA83D0BB, 0xBBA8E880941899E2), // -3.280306310840657572397879931169e-5
            F64x2::from_bits(0xBF0154ED4598D2CC, 0x3B825347BEB36A36), // -3.305768104803138394162995660892e-5
            F64x2::from_bits(0x3EA94F64F476E953, 0x3B4CADCA7FD7A90D), // 7.543007791079819166426271863463e-7
            F64x2::from_bits(0x3E9E127258537A6B, 0x3B2675D3546BB035), // 4.481085727601967420432483955913e-7
            F64x2::from_bits(0xBE4588C758E152E3, 0xBAC63124C21ABA65), // -1.002768603930923058466126457023e-8
        ],
        k_lo: [
            f64::from_bits(0xBE31021CDD89ACCF), // -3.96004248195682e-9
            f64::from_bits(0x3DD7CFA7A7922A26), // 8.662446618757066e-11
            f64::from_bits(0x3DBB0E010D0D1E73), // 2.4606109524669533e-11
            f64::from_bits(0xBD6276CC50CC5147), // -5.247800092124551e-13
            f64::from_bits(0xBD3FE925F3C66B00), // -1.1336970599743444e-13
            f64::from_bits(0x3CE5313E62607073), // 2.3528243973815007e-15
            f64::from_bits(0x3CBD04DB4AB85C22), // 4.027191237220428e-16
            f64::from_bits(0xBC62C068630536CC), // -8.132207067527432e-18
            f64::from_bits(0xBC34685E15C9731F), // -1.1063028544890084e-18
            f64::from_bits(0x3BDAB0CF4D8E0072), // 2.2607871107654313e-20
        ],
    },
    // [22.776546738526, 25.918139392115794)
    ZeroPoly {
        c: f64::from_bits(0x403854FA303820CA), // 2.433194257135691e1
        k_hi: [
            F64x2::from_bits(0x3C8B67DBFD3BD179, 0xB923BC27684365D3), // 4.754131206539137653892405856404e-17
            F64x2::from_bits(0xBFC4B2A38F1AB9B4, 0xBC63C998C5E747A1), // -1.617016266586240954888778571022e-1
            F64x2::from_bits(0x3F6B3878AADEB317, 0xBC06F6C520EF41C6), // 3.322826079019579430766873260167e-3
            F64x2::from_bits(0x3F9B750D89A9B35F, 0xBC22A69C6DCFB8CB), // 2.681370880820865434352740732922e-2
            F64x2::from_bits(0xBF41F6911725A944, 0x3B7CB257B02C6076), // -5.481918763455069679523360059528e-4
            F64x2::from_bits(0xBF55BEEE6FD51C8A, 0xBBEC4BAD037AAA9F), // -1.327259872896199284553944705433e-3
            F64x2::from_bits(0x3EFC3625D7A6507D, 0xBB83E8D7EA94ACD1), // 2.690459721461966435268168115456e-5
            F64x2::from_bits(0x3F005375A588A707, 0xBBA9BEC80AFE4282), // 3.113940028172221571240104456270e-5
            F64x2::from_bits(0xBEA4EE5E4E7CB465, 0x3B312D78519CF2BA), // -6.237961800208459701292102838030e-7
            F64x2::from_bits(0xBE9C7B3D81B5D0B9, 0x3AF344CFB7E81AB5), // -4.244060407201076613349961421594e-7
            F64x2::from_bits(0x3E41FCE14F49A01C, 0x3AE43926334EC40F), // 8.376228135753542414621365838035e-9
        ],
        k_lo: [
            f64::from_bits(0x3E30346643B86D26), // 3.7729473395589e-9
            f64::from_bits(0xBDD41C861C71C3B8), // -7.316492473601937e-11
            f64::from_bits(0xBDB9EEDA117C9850), // -2.3585939705575687e-11
            f64::from_bits(0x3D5F8CF2E5EE228A), // 4.483607274872538e-13
            f64::from_bits(0x3D3EC465294FC083), // 1.0930694074337473e-13
            f64::from_bits(0xBCE2501551059617), // -2.0331320253860417e-15
            f64::from_bits(0xBCBC21B6449E3636), // -3.9040558896199737e-16
            f64::from_bits(0x3C60612971E14493), // 7.10349255464326e-18
            f64::from_bits(0x3C33E0E3F645AFFE), // 1.077614474654969e-18
            f64::from_bits(0xBBD782E797A51D1B), // -1.9914884966514638e-20
        ],
    },
    // [25.918139392115794, 29.059732045705587)
    ZeroPoly {
        c: f64::from_bits(0x403B79ACEE8CFB7D), // 2.7475294980449224e1
        k_hi: [
            F64x2::from_bits(0x3C919E290514C619, 0xB92D5FEC5D3DE055), // 6.112406904655437179904531925605e-17
            F64x2::from_bits(0x3FC37AACEAC987B9, 0x3C2449CECCFBE447), // 1.521812578603752025233842990257e-1
            F64x2::from_bits(0xBF66AFE4FE0BC13E, 0x3C09E153C3470F03), // -2.769419909206898356745940954802e-3
            F64x2::from_bits(0xBF99DE7A33BC3A97, 0xBC2FDCBD2F73E546), // -2.526274624602339608268858295356e-2
            F64x2::from_bits(0x3F3E024F567AC4B6, 0xBBBB0C9B7DE6F12B), // 4.579013534381047651583734370023e-4
            F64x2::from_bits(0x3F548843C426ABDF, 0x3BE8512086CC2BF9), // 1.253191163963711209685834780398e-3
            F64x2::from_bits(0xBEF7A8E14711C114, 0x3B960D2799D3B61D), // -2.256363626823669891592799554770e-5
            F64x2::from_bits(0xBEFEECEB341AD7EB, 0x3B8AA756C7873806), // -2.949282066452139669811335873139e-5
            F64x2::from_bits(0x3EA1A743E05B3C55, 0x3B48B7A54D8FE9BD), // 5.261117059454345343965487509500e-7
            F64x2::from_bits(0x3E9B143D39C8DFE5, 0xBB3A797D057F5D7D), // 4.035094265694517321636354505716e-7
            F64x2::from_bits(0xBE3E8E000125623D, 0xBACD587D6C476791), // -7.114067573172209458751079097896e-9
        ],
        k_lo: [
            f64::from_bits(0xBE2EF28E321C527A), // -3.602761104106133e-9
            f64::from_bits(0x3DD137A1BFF346A1), // 6.263621564782915e-11
            f64::from_bits(0x3DB8E087B84DE702), // 2.26255641991121e-11
            f64::from_bits(0xBD5B3DCA3BC73696), // -3.871231099374144e-13
            f64::from_bits(0xBD3DA6A5499D8C98), // -1.0534136971790766e-13
            f64::from_bits(0x3CDFE6706C47A71D), // 1.7708142140615335e-15
            f64::from_bits(0x3CBB3B8E9A93D84A), // 3.77928864752977e-16
            f64::from_bits(0xBC5CC75C7A508D06), // -6.2403977057676726e-18
            f64::from_bits(0xBC33526F29D78929), // -1.0474481915265248e-18
            f64::from_bits(0x3BD4CD4976B7DB1D), // 1.7619896544002424e-20
        ],
    },
    // [29.059732045705587, 32.20132469929538)
    ZeroPoly {
        c: f64::from_bits(0x403E9E480605283C), // 3.0618286491641115e1
        k_hi: [
            F64x2::from_bits(0xBC719366994A86BC, 0xB8F250ED82394918), // -1.524456280251315087137212184560e-17
            F64x2::from_bits(0xBFC2740819F1CAAA, 0x3C6C99150772CD61), // -1.441660048181650375697273427001e-1
            F64x2::from_bits(0x3F6349369DC780CD, 0xBC0F3298E410CFEC), // 2.354246780882451748017829170123e-3
            F64x2::from_bits(0x3F98868D7401BF2E, 0x3C15B573B41488C1), // 2.395077724980792412412529152986e-2
            F64x2::from_bits(0xBF398CD1BEBE1451, 0x3BD9E57A62097EFE), // -3.898632119148692069207537570059e-4
            F64x2::from_bits(0xBF537EEF9AADEEE2, 0x3BF1A592993B1133), // -1.189931858532604732360635083188e-3
            F64x2::from_bits(0x3EF43394C95B2D36, 0xBB6A42D4ABBED026), // 1.926564126257282004538870404884e-5
            F64x2::from_bits(0x3EFD6DFCDB025FE3, 0xBB9B8A27FEF80653), // 2.806629135018713381935995964382e-5
            F64x2::from_bits(0xBE9E448FBC8A23A4, 0x3B33C61DE33DC388), // -4.510256386499434405506445363383e-7
            F64x2::from_bits(0xBE99D764EE079B52, 0x3B3E9E4706F5703E), // -3.850666257229861967293443526061e-7
            F64x2::from_bits(0x3E3A53958C8F6B8A, 0x3AAFF97F0260EEDA), // 6.129616100251015631164552620809e-9
        ],
        k_lo: [
            f64::from_bits(0x3E2DA0E1C87F32D7), // 3.4492049806383145e-9
            f64::from_bits(0xBDCDD7F919E64411), // -5.428527353106917e-11
            f64::from_bits(0xBDB7E68028B4937A), // -2.1737280850383682e-11
            f64::from_bits(0x3D57C2E68E01D5FB), // 3.376688153731973e-13
            f64::from_bits(0x3D3C97319246A8E0), // 1.0157421623586748e-13
            f64::from_bits(0xBCDC0354F5124465), // -1.5550347174193777e-15
            f64::from_bits(0xBCBA5AB6BDC6309D), // -3.657400897559826e-16
            f64::from_bits(0x3C5971EC78737DBD), // 5.5175080007599705e-18
            f64::from_bits(0x3C32C352F834D4D4), // 1.0171434422238744e-18
            f64::from_bits(0xBBD28030BD63741D), // -1.5670739471339238e-20
        ],
    },
    // [32.20132469929538, 35.34291735288517)
    ZeroPoly {
        c: f64::from_bits(0x4040E16907F8FB56), // 3.376101779610933e1
        k_hi: [
            F64x2::from_bits(0x3CBBEC2345AC18E2, 0xB95963D3D987B148), // 3.875013256161460457483896156443e-16
            F64x2::from_bits(0x3FC192F2627A74E3, 0xBC6E97BFA542387D), // 1.372969609118746692751899959458e-1
            F64x2::from_bits(0xBF60A846A83FEEAF, 0xBBFE500B58F61CA4), // -2.033365251916517348764675944042e-3
            F64x2::from_bits(0xBF975ECEAABF7F85, 0xBBE682B196CED8E0), // -2.282259862101867005018488569799e-2
            F64x2::from_bits(0x3F3617C581BE36D9, 0xBBD66BA20FFCD4CE), // 3.371102518422270695160730973940e-4
            F64x2::from_bits(0x3F529934B7A84481, 0x3BFE0AC8B6D39378), // 1.135159952017767189806674184280e-3
            F64x2::from_bits(0xBEF18123E875192B, 0xBB8D61FA5FCFEB21), // -1.669354835811505474843478180144e-5
            F64x2::from_bits(0xBEFC1F05A2D8511E, 0xBB531CDD3157D2B0), // -2.681844687486930929287122129443e-5
            F64x2::from_bits(0x3E9A4E0BC09269CB, 0xBB24D4F57848BCDB), // 3.919730607399489866916442049782e-7
            F64x2::from_bits(0x3E98BE81AD44AD5D, 0xBB39ACECDEBC155E), // 3.687168093231732245019694051082e-7
            F64x2::from_bits(0xBE36F73795E4173F, 0xBAD23AA262B7CDEC), // -5.347116831130590056372563238954e-9
        ],
        k_lo: [
            f64::from_bits(0xBE2C70AB156C4681), // -3.310864620124579e-9
            f64::from_bits(0x3DCA2621912F7189), // 4.7564662406801026e-11
            f64::from_bits(0x3DB701129BF33E77), // 2.092218910871245e-11
            f64::from_bits(0xBD54EC1687E27F1A), // -2.9732261158303776e-13
            f64::from_bits(0xBD3B990D2701CBDC), // -9.804728385367837e-14
            f64::from_bits(0x3CD8CBC581587941), // 1.3764535312123996e-15
            f64::from_bits(0x3CB983A9918FE011), // 3.540821269732378e-16
            f64::from_bits(0xBC56A54FFBC332BD), // -4.910514638260701e-18
            f64::from_bits(0xBC323776828CA913), // -9.875266874020211e-19
            f64::from_bits(0x3BD08B9AEF47EAD6), // 1.4014442060596843e-20
        ],
    },
    // [35.34291735288517, 38.48451000647496)
    ZeroPoly {
        c: f64::from_bits(0x404273A7B35A7AFF), // 3.690355531614295e1
        k_hi: [
            F64x2::from_bits(0x3CA334EB74FDCFD1, 0x394784BB69CFF9C8), // 1.332733769778246779282411711764e-16
            F64x2::from_bits(0xBFC0CF3EE98F769B, 0x3C50D95401A3F4C9), // -1.313246384078653147366908316611e-1
            F64x2::from_bits(0x3F5D26E7AF251E46, 0xBBD2054E6F67D637), // 1.779295209944368914152552171199e-3
            F64x2::from_bits(0x3F965D05948A946A, 0xBC09B55696A4EDE1), // 2.183922499775229463590949676173e-2
            F64x2::from_bits(0xBF335959B8482DDA, 0xBBDF04523FE57F77), // -2.952426951032518156417555478748e-4
            F64x2::from_bits(0xBF51CFF175D05C2B, 0x3BFB907FA5DB0178), // -1.087175179469000234974675162738e-3
            F64x2::from_bits(0x3EEEB594168790A2, 0xBB357172C6A15379), // 1.464333100015979513709042948588e-5
            F64x2::from_bits(0x3EFAF7544EEAC724, 0x3B8FD65D32FFB1B5), // 2.571690577707961208580373256647e-5
            F64x2::from_bits(0xBE9720522BB1FFAD, 0x3B3DCED12FB50AB2), // -3.446080360583425817233417265864e-7
            F64x2::from_bits(0xBE97C4126170295C, 0x3B1A0D132E4924C0), // -3.541395882661389719400389050380e-7
            F64x2::from_bits(0x3E343FA0EA5F89F1, 0xBAD532D5888B39BF), // 4.714482726059212440137350803864e-9
        ],
        k_lo: [
            f64::from_bits(0x3E2B5E23ABCC1C1C), // 3.1860233048832115e-9
            f64::from_bits(0xBDC722397F43628D), // -4.2079936677429457e-11
            f64::from_bits(0xBDB62F216683D904), // -2.0176324506256622e-11
            f64::from_bits(0x3D5294AB4644E008), // 2.640481744903562e-13
            f64::from_bits(0x3D3AACE12B7EE120), // 9.476974163675393e-14
            f64::from_bits(0xBCD61CDEF8371368), // -1.2275057218428872e-15
            f64::from_bits(0xBCB8B879448C4AA3), // -3.4306724682265614e-16
            f64::from_bits(0x3C54489FA23FE11B), // 4.398323246039918e-18
            f64::from_bits(0x3C31B10F91B9BA36), // 9.590659331023479e-19
            f64::from_bits(0xBBCDC1E6B0D80AD7), // -1.2602760255964043e-20
        ],
    },
    // [38.48451000647496, 41.62610266006476)
    ZeroPoly {
        c: f64::from_bits(0x404405E18393AFB5), // 4.0045944640266875e1
        k_hi: [
            F64x2::from_bits(0xBCAACA291AEE3C35, 0xB91831A6EA0719AB), // -1.858908146547320755058985689070e-16
            F64x2::from_bits(0x3FC0230BA90F2871, 0x3C6BCA82E12FBA6D), // 1.260695052608898201816718858958e-1
            F64x2::from_bits(0xBF59CA16F0C971A2, 0x3BF0C2667E91C092), // -1.574060824302760496750599477916e-3
            F64x2::from_bits(0xBF9579C1BDBCFC9A, 0x3C288939EC361D1F), // -2.097227783746955581786817643689e-2
            F64x2::from_bits(0x3F3120ECFAC5BF89, 0x3BB5BA72CA120054), // 2.613619388034329400455035951609e-4
            F64x2::from_bits(0x3F511DD26BBE2946, 0x3BF03A624D97B939), // 1.044707764750700508703658535608e-3
            F64x2::from_bits(0xBEEB37E7C9A570B9, 0xBB8B94C39391ACB9), // -1.297873523203211581520671227651e-5
            F64x2::from_bits(0xBEF9F01E7C190961, 0x3B9FB1CD5FBC48AC), // -2.473637119338714747040427888305e-5
            F64x2::from_bits(0x3E94887FE7A88EE5, 0xBB3C3D91CDD55FD3), // 3.059685479719845801667224165492e-7
            F64x2::from_bits(0x3E96E3723883F39E, 0x3AEBB166679861AC), // 3.410646560867967505732720577923e-7
            F64x2::from_bits(0xBE3204B644D800C0, 0x3AC788C42DBBE7AA), // -4.195237113865399785250936622095e-9
        ],
        k_lo: [
            f64::from_bits(0xBE2A659B13CFCF5D), // -3.073003322721631e-9
            f64::from_bits(0x3DC4A40C8818F1D1), // 3.7545425980185074e-11
            f64::from_bits(0x3DB56F0242DEF94E), // 1.9493771336237682e-11
            f64::from_bits(0xBD50A1015570E728), // -2.363112600032641e-13
            f64::from_bits(0xBD39D22289E5347F), // -9.173404975381921e-14
            f64::from_bits(0x3CD3DB238F7F2033), // 1.1022300494922903e-15
            f64::from_bits(0x3CB7F9CC81E0C782), // 3.327307427400872e-16
            f64::from_bits(0xBC52473C8445DCF6), // -3.963467393137084e-18
            f64::from_bits(0xBC313139747765F8), // -9.319955258346931e-19
            f64::from_bits(0x3BCAE75A21C0B31E), // 1.1394168220712376e-20
        ],
    },
    // [41.62610266006476, 44.76769531365455)
    ZeroPoly {
        c: f64::from_bits(0x4045981787D668DB), // 4.318821809739321e1
        k_hi: [
            F64x2::from_bits(0x3C8EF92AA9C8E53C, 0x391D1090C39E51F4), // 5.373012319630951672687347743303e-17
            F64x2::from_bits(0xBFBF13FB0C0E6FCD, 0x3C463CEA806B69AD), // -1.213986305651239664719352814681e-1
            F64x2::from_bits(0x3F5706ED3D935C85, 0xBBF62588C13A8B89), // 1.405460052685631909792791218794e-3
            F64x2::from_bits(0x3F94AF74CBD77BF0, 0xBC3DF2D80FF3BB9D), // 2.020056242133344397167777808197e-2
            F64x2::from_bits(0xBF2E9A9E66E5A740, 0x3B90C18A2E98BD30), // -2.334898339982297695860535334191e-4
            F64x2::from_bits(0xBF507EC9ED824FCB, 0x3BF6C4A7B9454C59), // -1.006791300857112302307734059429e-3
            F64x2::from_bits(0x3EE856D4518AB279, 0xBB62DFDBF75E82F0), // 1.160582410028001084528433582272e-5
            F64x2::from_bits(0x3EF9040DE830645F, 0xBB921F0E80083C66), // 2.385696144466644025616880735650e-5
            F64x2::from_bits(0xBE9262F69C56C8FB, 0xBB289FE9FB8BBA77), // -2.739813250398663680274781206677e-7
            F64x2::from_bits(0xBE9618C94A544AA7, 0x3AFE149D5AED4C6A), // -3.292682982529840730515728262903e-7
            F64x2::from_bits(0x3E3029D2C8C0AEFB, 0xBAB663AFE2A08A93), // 3.763328437110108127397126404199e-9
        ],
        k_lo: [
            f64::from_bits(0x3E2983BCA085BAE7), // 2.9702900175730903e-9
            f64::from_bits(0xBDC28E2A00A6DEBF), // -3.375194575345285e-11
            f64::from_bits(0xBDB4BEEABD4FF0B6), // -1.886816729941623e-11
            f64::from_bits(0x3D4DF8814B1CB12E), // 2.1295479414233272e-13
            f64::from_bits(0x3D3907B691EA3D29), // 8.892488362313906e-14
            f64::from_bits(0xBCD1F12262AD1B5D), // -9.959772412575452e-16
            f64::from_bits(0xBCB7477E78B6E353), // -3.23064818682511e-16
            f64::from_bits(0x3C509075D37AAAD2), // 3.5918095501891046e-18
            f64::from_bits(0x3C30B85FD61B5C95), // 9.064045275369398e-19
            f64::from_bits(0xBBC8729AC64DB612), // -1.0353992754273785e-20
        ],
    },
    // [44.76769531365455, 47.909287967244346)
    ZeroPoly {
        c: f64::from_bits(0x40472A4A85CC317E), // 4.633039925070169e1
        k_hi: [
            F64x2::from_bits(0x3CB13E3C695EA088, 0x393B061DCBC3A2F4), // 2.392962121214805367699240482680e-16
            F64x2::from_bits(0x3FBE018DAC1C17E3, 0xBC5EB9933038F4A0), // 1.172112031163972462033868123234e-1
            F64x2::from_bits(0xBF54B994DD05C422, 0x3BDEB369470CABCA), // -1.264949201950090326376373854516e-3
            F64x2::from_bits(0xBF93F9E0DB07E7EF, 0x3C285EC817AEC59A), // -1.950789772852340040468035593173e-2
            F64x2::from_bits(0x3F2B8E55B75B1519, 0x3BC418D4FB765680), // 2.102355608057181113523354400079e-4
            F64x2::from_bits(0x3F4FDF68A78BB3D0, 0x3BE2CDFB1787C5AA), // 9.726773269825798497940213756252e-4
            F64x2::from_bits(0xBEE5EE9D17106AD1, 0x3B81A8D532AE7769), // -1.045803284865001707443892437205e-5
            F64x2::from_bits(0xBEF82EE6DFDFEDAC, 0xBB7FE8474874CA21), // -2.306290661378627319412088578883e-5
            F64x2::from_bits(0x3E90962D7F6D6495, 0x3B2D3C30F2C92046), // 2.471600732321697290068005767950e-7
            F64x2::from_bits(0x3E9560EDCE7D4082, 0x3B36637D2CE8F441), // 3.185663914788832879619706109400e-7
            F64x2::from_bits(0xBE2D34381D050018, 0x3AC19F881C918E9F), // -3.399790872573670464503028585319e-9
        ],
        k_lo: [
            f64::from_bits(0xBE28B5A0E72D71B4), // -2.8765628155674476e-9
            f64::from_bits(0x3DC0CAB95487A13E), // 3.054427072919772e-11
            f64::from_bits(0x3DB41D21B936F464), // 1.8293390739648633e-11
            f64::from_bits(0xBD4B2D2E35CCD86F), // -1.9310054967318525e-13
            f64::from_bits(0xBD384C4BE9436C49), // -8.632395531819291e-14
            f64::from_bits(0x3CD04E285D6152D5), // 9.051261638778547e-16
            f64::from_bits(0x3CB6A0FF6299EB95), // 3.140390290632327e-16
            f64::from_bits(0xBC4E2DD5EAC9BE48), // -3.2720186549929565e-18
            f64::from_bits(0xBC3046857A8E68FA), // -8.822952255875517e-19
            f64::from_bits(0x3BC6526F0A5809C8), // 9.453737622060176e-21
        ],
    },
];

#[cfg(test)]
mod tests {
    use super::{
        AIRY_T_MIN, AIRY_TBL, DEBYE_A1, DEBYE_B1, DEBYE_BELOW_MIN, DEBYE_DIST_MIN, DEBYE_FAR_R_MIN,
        DEBYE_POLY, DEBYE_TERMS, F64x2, FRAC_1_PI, J0_LARGE_MIN, J0_SMALL_MAX, J0_TBL,
        J1_LARGE_MIN, J1_SMALL_MAX, J1_TBL, LARGE0, LARGE1, Large, ORDER_DEBYE_FAR_MIN,
        ORDER_DEBYE_MIN, ORDER_LARGE_FACTOR, PHASE_FIXED_BITS, Y0_ANCHORED, Y0_LARGE_MIN,
        Y0_SMALL_MAX, Y0_TBL, Y1_ANCHORED, Y1_LARGE_MIN, Y1_SMALL_MAX, Y1_TBL, ZeroPoly, airy,
        atan_precise, atanh_precise, debye_above_parts, debye_below, debye_far, debye_valid,
        jn_backward_terms, jn_debye, jn_log2_bound, jn_underflows, olver, order_large_parts,
        phase_precise, yn_debye, yn_overflows,
    };
    use crate::traits::Float as _;
    use rug::ops::Pow as _;

    const PREC: u32 = 512;

    fn float(x: f64) -> rug::Float {
        rug::Float::with_val(PREC, x)
    }

    fn pi() -> rug::Float {
        rug::Float::with_val(PREC, rug::float::Constant::Pi)
    }

    /// Returns `log2` of the relative error of `approx` with respect to
    /// `exact`.
    fn rel_error_log2(approx: &rug::Float, exact: &rug::Float) -> f64 {
        let err = rug::Float::with_val(PREC, approx - exact) / exact;
        err.abs().log2().to_f64()
    }

    /// A table of polynomials and how its intervals are selected.
    struct Table {
        tbl: &'static [ZeroPoly],
        /// The functions that the polynomials approximate
        func: fn(&rug::Float) -> rug::Float,
        /// The first `anchored` polynomials approximate
        /// `scale(x) * (func(x) - 2/π * first_kind(x) * ln(x / c))`
        anchored: usize,
        first_kind: fn(&rug::Float) -> rug::Float,
        scale: fn(&rug::Float) -> rug::Float,
        /// The index of the interval of `x` is
        /// `floor(x / π + index_offset) - 1`
        index_offset: f64,
        /// The table is used for `x_min <= x < x_max`
        x_min: f64,
        x_max: f64,
    }

    fn tables() -> [Table; 4] {
        let one = |_: &rug::Float| float(1.0);
        let ident = |x: &rug::Float| x.clone();
        [
            Table {
                tbl: &J0_TBL,
                func: |x| x.clone().j0(),
                anchored: 0,
                first_kind: |x| x.clone().j0(),
                scale: one,
                index_offset: 0.75,
                x_min: J0_SMALL_MAX,
                x_max: J0_LARGE_MIN,
            },
            Table {
                tbl: &J1_TBL,
                func: |x| x.clone().j1(),
                anchored: 0,
                first_kind: |x| x.clone().j1(),
                scale: one,
                index_offset: 0.25,
                x_min: J1_SMALL_MAX,
                x_max: J1_LARGE_MIN,
            },
            Table {
                tbl: &Y0_TBL,
                func: |x| x.clone().y0(),
                anchored: Y0_ANCHORED,
                first_kind: |x| x.clone().j0(),
                scale: one,
                index_offset: 1.25,
                x_min: Y0_SMALL_MAX,
                x_max: Y0_LARGE_MIN,
            },
            Table {
                tbl: &Y1_TBL,
                func: |x| x.clone().y1(),
                anchored: Y1_ANCHORED,
                first_kind: |x| x.clone().j1(),
                scale: ident,
                index_offset: 0.75,
                x_min: Y1_SMALL_MAX,
                x_max: Y1_LARGE_MIN,
            },
        ]
    }

    #[test]
    fn test_zero_polys() {
        // `ZeroPoly::eval` requires a relative error less than 2^-70 from the
        // polynomials, terms of degree 11 and higher that are less than
        // 2^-16 relative to the result, and exact differences with the
        // centers, which must be very close to the zeros.
        let frac_2_pi = 2 / pi();
        for (t, table) in tables().iter().enumerate() {
            // The limits of the table select its first and last intervals
            let index = |x: f64| (x * FRAC_1_PI + table.index_offset) as usize - 1;
            let below_max = f64::from_bits(table.x_max.to_bits() - 1);
            assert_eq!(index(table.x_min), 0);
            assert_eq!(index(below_max), table.tbl.len() - 1);

            for (i, poly) in table.tbl.iter().enumerate() {
                // The interval is `[a, b)`
                let limit = |k: usize| (float(k as f64 + 1.0 - table.index_offset) * pi()).to_f64();
                let a = limit(i).max(table.x_min);
                let b = limit(i + 1).min(table.x_max);
                let c = float(poly.c);
                assert!(a < poly.c && poly.c < b);

                // The center is next to the zero, so `K0` is tiny
                let k1 = poly.k_hi[1].hi();
                assert!(
                    poly.k_hi[0].hi().abs() <= (2.0 * f64::EPSILON * poly.c) * k1.abs(),
                    "table = {t}, i = {i}"
                );

                // The function that the polynomial approximates, and its
                // contribution to the result
                let anchored = i < table.anchored;
                let exact = |x: &rug::Float| {
                    let f = (table.func)(x);
                    if anchored {
                        let ln = rug::Float::with_val(PREC, x / &c).ln();
                        (f - ln * (table.first_kind)(x) * &frac_2_pi) * (table.scale)(x)
                    } else {
                        f
                    }
                };
                let result = |x: &rug::Float| {
                    if anchored {
                        (table.func)(x) * (table.scale)(x)
                    } else {
                        (table.func)(x)
                    }
                };

                // Values of the interval, and some beyond its limits, where
                // the rounding of the index can select this interval too
                let width = b - a;
                let xs = (0..=128)
                    .map(|j| a + width * (f64::from(j) / 128.0))
                    .chain([a - 1e-12 * a, b + 1e-12 * b, poly.c])
                    // The values closest to the center, whose results are
                    // the smallest ones
                    .chain((1..=4).flat_map(|j| {
                        let bits = poly.c.to_bits();
                        [f64::from_bits(bits - j), f64::from_bits(bits + j)]
                    }))
                    // Without excess precision (x87)
                    .map(f64::purify)
                    .filter(|&x| x >= table.x_min && x < table.x_max);
                for x in xs {
                    // `d` is exact
                    let d = x - poly.c;
                    let xr = float(x);
                    assert_eq!(float(d), rug::Float::with_val(PREC, &xr - &c), "x = {x:e}");

                    let coeffs = (poly.k_hi.iter().map(|k| float(k.hi()) + k.lo()))
                        .chain(poly.k_lo.iter().map(|&k| float(k)));
                    let approx = coeffs
                        .clone()
                        .rev()
                        .fold(rug::Float::new(PREC), |acc, k| acc * d + k);
                    let res = result(&xr);
                    let err = rug::Float::with_val(PREC, &approx - exact(&xr)) / &res;
                    if !err.is_zero() {
                        let err = err.abs().log2().to_f64();
                        assert!(
                            err < -70.0,
                            "table = {t}, i = {i}, x = {x:e}, err = 2^{err}"
                        );
                    }

                    let tail = coeffs
                        .rev()
                        .take(poly.k_lo.len())
                        .fold(rug::Float::new(PREC), |acc, k| acc * d + k)
                        * float(d).pow(poly.k_hi.len() as u32);
                    let ratio = (tail / res).abs();
                    let ratio = if ratio.is_zero() {
                        f64::NEG_INFINITY
                    } else {
                        ratio.log2().to_f64()
                    };
                    assert!(ratio < -16.0, "table = {t}, i = {i}, x = {x:e}");
                }
            }
        }
    }

    /// Returns `(modulus, phase)` of the Bessel functions of order `nu`,
    /// where the phase is the term that tends to zero:
    /// `θ(x) - x + (2 * nu + 1) * π/4`.
    fn modulus_phase(nu: i32, x: &rug::Float) -> (rug::Float, rug::Float) {
        let j = x.clone().jn(nu);
        let y = x.clone().yn(nu);
        // (j + i * y) * exp(-i * (x - (2 * nu + 1) * π/4)) = M * exp(i * phase)
        let w: rug::Float = x.clone() - pi() * (2 * nu + 1) / 4;
        let (s, c) = w.sin_cos(rug::Float::new(PREC));
        let re = rug::Float::with_val(PREC, &j * &c) + rug::Float::with_val(PREC, &y * &s);
        let im = rug::Float::with_val(PREC, &y * &c) - rug::Float::with_val(PREC, &j * &s);
        (re.clone().hypot(&im), im.atan2(&re))
    }

    fn larges() -> [(i32, &'static Large); 2] {
        [(0, &LARGE0), (1, &LARGE1)]
    }

    #[test]
    fn test_large_polys() {
        // `Large::eval` requires a relative error less than 2^-78 from the
        // polynomial of the modulus and an absolute error less than 2^-82.5
        // from the one of the phase
        for (nu, large) in larges() {
            assert_eq!(i32::from(large.offset), (2 * nu + 1) * 32);

            let xs = (0..=64)
                .map(|j| J0_LARGE_MIN * (1.0 + f64::from(j) / 64.0))
                .chain((7..35).map(|e| f64::exp2i_fast(e) * 1.25));
            for x in xs {
                let xr = float(x);
                let s = rug::Float::with_val(PREC, xr.square_ref()).recip();
                let (modulus, phase) = modulus_phase(nu, &xr);

                let (a1, a) = large.modulus;
                let poly = a
                    .iter()
                    .rev()
                    .fold(rug::Float::new(PREC), |acc, &k| acc * &s + k);
                let poly = (poly * &s + a1.hi() + a1.lo()) * &s + 1;
                let scale: rug::Float = 2 / (pi() * &xr);
                let approx = scale.sqrt() * poly;
                assert!(rel_error_log2(&approx, &modulus) < -78.0, "x = {x:e}");

                let ([b0, b1], b) = large.phase;
                let poly = b
                    .iter()
                    .rev()
                    .fold(rug::Float::new(PREC), |acc, &k| acc * &s + k);
                let poly = (poly * &s + b1.hi() + b1.lo()) * &s + b0.hi() + b0.lo();
                let err: rug::Float = poly / &xr - phase;
                assert!(err.abs().log2().to_f64() < -82.5, "x = {x:e}");
            }
        }
    }

    #[test]
    fn test_phase_precise() {
        // The phase is `n * π/128 + b`, with an absolute error less than
        // 2^-128, or 2^-100 relative to `b` when it is not small (as a sum of
        // two `f64`)
        for (nu, large) in larges() {
            assert_eq!(PHASE_FIXED_BITS, 123);

            // Arguments from the smallest one, with all the exponents
            let spread = (0..=512)
                .map(|j| J0_LARGE_MIN * (1.0 + f64::from(j) / 512.0))
                .chain((6..1024).flat_map(|e| {
                    [1.0, 1.2345678901234567, 1.9999999999999998].map(|m| f64::exp2i_fast(e) * m)
                }));
            // The arguments closest to some zeros of the functions of both
            // kinds, where `b` is tiny
            let zeros = (6..120).flat_map(|e| {
                let prec = PREC + 256;
                let pi = rug::Float::with_val(prec, rug::float::Constant::Pi);
                [0, 1].map(|second| {
                    // The zeros are close to `(k + (2 * nu - 1 - 2 * second) / 4) * π`
                    let x0 = rug::Float::with_val(prec, f64::exp2i_fast(e) * 1.37);
                    let offset = rug::Float::with_val(prec, 2 * nu - 1 - 2 * second) / 4u32;
                    let k = (x0 / &pi - &offset).round();
                    let mut z = (k + offset) * &pi;
                    for _ in 0..8 {
                        let (f, f_next) = if second == 0 {
                            (z.clone().jn(nu), z.clone().jn(nu + 1))
                        } else {
                            (z.clone().yn(nu), z.clone().yn(nu + 1))
                        };
                        // f'(z) = nu / z * f(z) - f_next(z)
                        let df = rug::Float::with_val(prec, &f * nu) / &z - f_next;
                        z -= f / df;
                    }
                    z.to_f64()
                })
            });
            for x in spread.chain(zeros) {
                let (n, b) = phase_precise(x, large);
                assert!(b.bh.abs() <= 0.0123, "x = {x:e}");

                // The precision must be enough for the reduction of the
                // largest arguments
                let prec = PREC + 1100;
                let xr = rug::Float::with_val(prec, x);
                let pi = rug::Float::with_val(prec, rug::float::Constant::Pi);
                let j = xr.clone().jn(nu);
                let y = xr.clone().yn(nu);
                // θ(x) = atan2(y, j) (modulo 2 * π)
                let theta = y.atan2(&j);
                let approx: rug::Float =
                    rug::Float::with_val(prec, u32::from(n)) * &pi / 128u32 + b.bh + b.bl;
                // The difference is a multiple of `2 * π`
                let diff: rug::Float = (approx - theta) / (2u32 * pi.clone());
                let frac: rug::Float = diff.clone() - diff.round();
                let err = (frac * 2u32 * pi).abs().to_f64();
                let limit = (b.bh.abs() * f64::exp2i_fast(-100)).max(f64::exp2i_fast(-128));
                assert!(err < limit, "x = {x:e}, b = {:e}, err = {err:e}", b.bh);
            }
        }
    }

    #[test]
    fn test_order_large_parts() {
        // The modulus has a relative error less than 2^-74, and the phase an
        // absolute error less than 2^-87. The functions of `rug` are not
        // accurate from the order 2^30 with 32-bit targets, and the largest
        // order here is the one whose square is not exact as `f64`.
        let orders: [i32; 8] = [2, 3, 4, 7, 40, 1000, 65_536, (1 << 29) + 3];
        for n in orders {
            // The smallest arguments, where the truncation errors are the
            // largest, and arguments with all the larger exponents
            let nf = f64::from(n);
            let x_min = nf * nf * ORDER_LARGE_FACTOR;
            let factors = [1.0, 1.0 + f64::EPSILON, 1.37, 2.0 - f64::EPSILON];
            let xs = factors
                .map(|f| x_min * f)
                .into_iter()
                .chain((x_min.exponent() + 1..1024).map(|e| f64::exp2i_fast(e) * 1.37));
            for x in xs {
                let (m_hi, m_lo, k, b) = order_large_parts(n.unsigned_abs(), x);
                assert!(b.bh.abs() <= 0.0123, "n = {n}, x = {x:e}");

                // The precision must be enough for the reduction of the
                // largest arguments
                let prec = PREC + 1100;
                let xr = rug::Float::with_val(prec, x);
                let pi = rug::Float::with_val(prec, rug::float::Constant::Pi);
                let j = xr.clone().jn(n);
                let y = xr.clone().yn(n);

                let modulus = j.clone().hypot(&y);
                let approx = rug::Float::with_val(prec, m_hi) + m_lo;
                let err = rel_error_log2(&approx, &modulus);
                assert!(err < -74.0, "n = {n}, x = {x:e}, modulus error = 2^{err}");

                // θ(x) = atan2(y, j) (modulo 2 * π)
                let theta = y.atan2(&j);
                let approx: rug::Float =
                    rug::Float::with_val(prec, u32::from(k)) * &pi / 128u32 + b.bh + b.bl;
                // The difference is a multiple of `2 * π`
                let diff: rug::Float = (approx - theta) / (2u32 * pi.clone());
                let frac: rug::Float = diff.clone() - diff.round();
                let err = (frac * 2u32 * pi).abs().log2().to_f64();
                assert!(err < -87.0, "n = {n}, x = {x:e}, phase error = 2^{err}");
            }
        }
    }

    #[test]
    fn test_atan_atanh_precise() {
        // Values around the points of the tables and between them (the
        // largest `u`), with and without a low part
        for hyperbolic in [false, true] {
            let max = if hyperbolic { 0.9375 } else { 1.0 };
            let ys = (1..=1024)
                .flat_map(|i| [0.0, 0.3, 0.5, 0.9999].map(|f| (f64::from(i) + f) / 1024.0 * max));
            for y in ys {
                for y_lo in [0.0, y * f64::exp2i_fast(-54) * 0.7] {
                    let v = F64x2::add11(y, y_lo);
                    let r = if hyperbolic {
                        atanh_precise(v)
                    } else {
                        atan_precise(v)
                    };
                    let yr = float(v.hi()) + v.lo();
                    let exact = if hyperbolic { yr.atanh() } else { yr.atan() };
                    let err = rel_error_log2(&(float(r.hi()) + r.lo()), &exact);
                    // The documented errors (see `atan_precise` and
                    // `atanh_precise`)
                    let limit = if hyperbolic && y > 0.75 { -89.0 } else { -99.0 };
                    assert!(
                        err < limit,
                        "y = {y:e}, hyperbolic = {hyperbolic}, err = 2^{err}"
                    );
                }
            }
        }
    }

    /// Returns `jn(n, x)` (`second = false`) or `yn(n, x)`, from `rug` (`jn`
    /// with `x < n`) or from the forward recurrence of the functions of orders
    /// 0 and 1 (which is stable with the other arguments), with 512 bits.
    fn order_rug(n: u32, x: f64, second: bool) -> rug::Float {
        let xr = float(x);
        if !second && x < f64::from(n) {
            return xr.jn(n as i32);
        }
        let (mut a, mut b) = if second {
            (xr.clone().y0(), xr.clone().y1())
        } else {
            (xr.clone().j0(), xr.clone().j1())
        };
        for k in 1..n {
            let c = rug::Float::with_val(PREC, &b * (2 * k)) / &xr - &a;
            a = b;
            b = c;
        }
        b
    }

    /// Returns the arguments `x < n` (`above = false`) or `x > n` with which
    /// the Debye expansions are used for the order `n`, from the closest ones
    /// to `n` (where their truncation errors are the largest). Below the
    /// order, they go down to the first ones where they are used (see
    /// `debye_below_used`).
    fn debye_args(n: u32, above: bool) -> impl Iterator<Item = f64> {
        let nf = f64::from(n);
        let n_r = float(nf);
        // R^3 / n^2 = DEBYE_DIST_MIN
        let r_min = rug::Float::with_val(PREC, n_r.square_ref()) * DEBYE_DIST_MIN;
        let r_min = r_min.cbrt();
        let factors = if above {
            [1.0, 1.0001, 1.05, 1.5, 3.0, 30.0, 1e3, 1e5]
        } else {
            [1.0, 1.0001, 1.05, 1.5, 3.0, 10.0, 30.0, f64::NAN]
        };
        let xs = factors.map(|f| {
            let r2 = rug::Float::with_val(PREC, (r_min.clone() * f).square_ref());
            let n2 = rug::Float::with_val(PREC, n_r.square_ref());
            if above { n2 + r2 } else { n2 - r2 }.sqrt().to_f64()
        });
        // Around `R = n` (the phase is calculated differently on each side),
        // and close to the limit of `order_large_eval`
        let s2 = core::f64::consts::SQRT_2 * nf;
        let extra = if above {
            [
                s2 * (1.0 - 1e-9),
                s2 * (1.0 + 1e-9),
                nf * nf * ORDER_LARGE_FACTOR * 0.999,
            ]
        } else {
            // Where `debye_below` starts being used
            let first = |pred: &dyn Fn(f64) -> bool| {
                let (mut lo, mut hi) = (0.25 * nf, nf);
                while f64::from_bits(lo.to_bits() + 1) < hi {
                    let mid = 0.5 * (lo + hi);
                    if pred(mid) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                hi
            };
            let x_min = (DEBYE_BELOW_MIN * nf).purify();
            [
                first(&|x| jn_underflows(n, x, 1075)).max(x_min),
                first(&|x| yn_overflows(n, x, 1024)).max(x_min),
                x_min,
            ]
        };
        xs.into_iter()
            .chain(extra)
            .map(f64::purify)
            .filter(move |&x| {
                debye_valid(nf, x)
                    && (above || debye_below_used(n, x, false) || debye_below_used(n, x, true))
            })
    }

    /// Returns `true` if `debye_below` is used with `x < n` for the functions
    /// of the first kind (`second = false`) or of the second one, which is
    /// when `x >= DEBYE_BELOW_MIN * n` and `jn_underflows` or `yn_overflows`
    /// is false.
    fn debye_below_used(n: u32, x: f64, second: bool) -> bool {
        let checked = if second {
            !yn_overflows(n, x, 1024)
        } else {
            !jn_underflows(n, x, 1075)
        };
        checked && x >= DEBYE_BELOW_MIN * f64::from(n)
    }

    #[test]
    fn test_debye() {
        // The results are compared before the final rounding, with limits
        // somewhat above the documented errors (without the one of the sine)
        for n in [ORDER_DEBYE_MIN, 2000, 4096, 20_011] {
            for x in debye_args(n, false) {
                for second in [false, true] {
                    if !debye_below_used(n, x, second) {
                        continue;
                    }
                    let e = debye_below(n, x, second);
                    let approx = (float(e.hi()) + e.lo()).exp();
                    let approx = if second { -approx } else { approx };
                    let exact = order_rug(n, x, second);
                    let err = rel_error_log2(&approx, &exact);
                    assert!(
                        err < -65.0,
                        "n = {n}, x = {x:e}, second = {second}, err = 2^{err}"
                    );
                }
            }
            for x in debye_args(n, true) {
                check_debye_above(n, x);
            }
        }
    }

    /// Checks the modulus and the phase of `debye_above_parts`, with limits
    /// somewhat above their documented errors.
    fn check_debye_above(n: u32, x: f64) {
        let (m_hi, m_lo, k, c) = debye_above_parts(n, x);
        let j = order_rug(n, x, false);
        let y = order_rug(n, x, true);
        let modulus = j.clone().hypot(&y);
        let err = rel_error_log2(&(float(m_hi) + m_lo), &modulus);
        assert!(err < -68.0, "n = {n}, x = {x:e}, modulus error = 2^{err}");

        // The phase error, from M * exp(i * θ) = j + i * y
        let theta = float(f64::from(k)) * pi() / 128u32 + c.bh + c.bl;
        let (s, co) = theta.sin_cos(rug::Float::new(PREC));
        let err_j = rug::Float::with_val(PREC, &modulus * &co) - &j;
        let err_y = rug::Float::with_val(PREC, &modulus * &s) - &y;
        let err = (err_j.hypot(&err_y) / &modulus).log2().to_f64();
        assert!(err < -66.0, "n = {n}, x = {x:e}, phase error = 2^{err}");
    }

    #[test]
    fn test_debye_far() {
        // The orders below `ORDER_DEBYE_MIN` with which `debye_above` is used
        // far beyond them, from the closest arguments, where the truncation
        // error is the largest, up to where `order_large_eval` is used
        for n in [ORDER_DEBYE_FAR_MIN, 21, 100, 600, ORDER_DEBYE_MIN - 1] {
            let nf = f64::from(n);
            // The first argument with which it is used, which the rounding
            // errors only move slightly from `sqrt(n^2 + DEBYE_FAR_R_MIN^2)`
            let x_limit = crate::generic::sqrt(nf * nf + DEBYE_FAR_R_MIN * DEBYE_FAR_R_MIN);
            let mut x_min = f64::from_bits(x_limit.to_bits() - 4);
            while !debye_far(n, x_min) {
                x_min = f64::from_bits(x_min.to_bits() + 1);
            }
            assert!(x_min.to_bits().abs_diff(x_limit.to_bits()) <= 2, "n = {n}");

            let x_max = nf * nf * ORDER_LARGE_FACTOR;
            let x_next = f64::from_bits(x_min.to_bits() + 1);
            let xs = [1.0, 1.0001, 1.05, 1.5, 3.0, 30.0, 1e3, 1e5].map(|f| x_min * f);
            for x in xs.into_iter().chain([x_next, x_max * 0.999]) {
                if x < x_max {
                    assert!(debye_far(n, x), "n = {n}, x = {x:e}");
                    check_debye_above(n, x);
                }
            }
        }
        let x = 1e4;
        assert!(debye_far(ORDER_DEBYE_FAR_MIN, x) && !debye_far(ORDER_DEBYE_FAR_MIN - 1, x));
        assert!(!debye_far(ORDER_DEBYE_MIN, 1e5));
    }

    /// Returns `(Ai(t), Ai'(t))` (`second = false`) or `(Bi(t), Bi'(t))`, from
    /// their Maclaurin series (`rug` only has `Ai`), for `|t| < 32`.
    ///
    /// Their terms are up to about 2^174 times larger than the smallest
    /// results (`Ai(31)` is about 2^-170), which the precision makes up for.
    fn airy_rug(t: &rug::Float, second: bool) -> (rug::Float, rug::Float) {
        let prec = 1024;
        let three = rug::Float::with_val(prec, 3);
        let third = rug::Float::with_val(prec, three.recip_ref());
        let two_thirds = rug::Float::with_val(prec, &third * 2u32);
        // Ai(0) = 3^(-2/3) / Γ(2/3) and Ai'(0) = -3^(-1/3) / Γ(1/3), and the
        // ones of `Bi` are `sqrt(3)` and `-sqrt(3)` times them
        let mut c0 = rug::Float::with_val(prec, (&three).pow(-two_thirds.clone()));
        c0 /= two_thirds.gamma();
        let mut c1 = -rug::Float::with_val(prec, (&three).pow(-third.clone()));
        c1 /= third.gamma();
        if second {
            let sqrt3 = three.sqrt();
            c0 *= &sqrt3;
            c1 *= -sqrt3;
        }

        // y = sum(c[m] * t^m) and y' = sum(m * c[m] * t^(m - 1)), with
        // (m + 2) * (m + 1) * c[m + 2] = c[m - 1]. The terms are less than
        // 2^-500 times the results from m = 1000.
        let mut y = c0.clone();
        let mut dy = rug::Float::new(prec);
        // c[m % 3] has c[m] from c[m - 3]
        let mut c = [c0, c1, rug::Float::new(prec)];
        // t^(m - 1)
        let mut t_pow = rug::Float::with_val(prec, 1);
        for m in 1..1000u32 {
            let cm = &mut c[(m % 3) as usize];
            if m >= 3 {
                *cm /= (m - 1) * m;
            }
            dy += rug::Float::with_val(prec, &t_pow * &*cm) * m;
            t_pow *= t;
            y += rug::Float::with_val(prec, &t_pow * &*cm);
        }
        (y, dy)
    }

    #[test]
    fn test_airy() {
        // `airy_rug`, checked with the `Ai` of `rug`. Its `Bi` is compared
        // below with `AIRY_TBL`, which the generator checks independently
        // (with the asymptotic expansions).
        for t in [-25.0, -1.5, 0.0, 2.0, 30.0] {
            let (ai, _) = airy_rug(&float(t), false);
            let err = rel_error_log2(&ai, &float(t).ai());
            assert!(err < -200.0, "t = {t}, err = 2^{err}");
        }

        // At the points of `AIRY_TBL`, between them and up to 1/8 beyond the
        // first and the last ones, also with the low part of `t`. The errors
        // are relative to the results with `t > 0`, and to the modulus of the
        // functions (or of their derivatives) otherwise.
        let step = 1.0 / 16.0;
        let t_max = AIRY_T_MIN + (AIRY_TBL.len() - 1) as f64 / 4.0;
        let count = ((t_max - AIRY_T_MIN + 0.25) / step) as i32;
        for i in 0..=count {
            let t_hi = AIRY_T_MIN - 0.125 + f64::from(i) * step;
            let t_lo = t_hi * f64::exp2i_fast(-60);
            let t = float(t_hi) + t_lo;
            let (ai, dai) = airy_rug(&t, false);
            let (bi, dbi) = airy_rug(&t, true);
            for second in [false, true] {
                let (y_exact, dy_exact) = if second { (&bi, &dbi) } else { (&ai, &dai) };
                let (scale, scale_d) = if t_hi > 0.0 {
                    (y_exact.clone().abs(), dy_exact.clone().abs())
                } else {
                    (ai.clone().hypot(&bi), dai.clone().hypot(&dbi))
                };
                let (y, dy) = airy(F64x2::new(t_hi, t_lo), second);
                let err = rug::Float::with_val(PREC, float(y.hi()) + y.lo() - y_exact) / &scale;
                let err = err.abs().log2().to_f64();
                let err_d = rug::Float::with_val(PREC, float(dy) - dy_exact) / &scale_d;
                let err_d = err_d.abs().log2().to_f64();
                assert!(err < -75.0, "t = {t_hi}, second = {second}, err = 2^{err}");
                assert!(
                    err_d < -48.0,
                    "t = {t_hi}, second = {second}, err_d = 2^{err_d}"
                );
            }
        }
    }

    #[test]
    fn test_olver_args() {
        // `t` and `ζ` of `olver` at the limits where `debye_valid` turns
        // false, where they are the largest in magnitude
        let (mut t_min, mut t_max) = (f64::INFINITY, f64::NEG_INFINITY);
        let (mut z_min, mut z_max) = (f64::INFINITY, f64::NEG_INFINITY);
        for n in [
            ORDER_DEBYE_MIN,
            ORDER_DEBYE_MIN + 1,
            1100,
            2000,
            4096,
            1 << 20,
            (1 << 31) - 1,
            1 << 31,
        ] {
            let nf = f64::from(n);
            for above in [false, true] {
                // The bisection keeps `lo` valid and `hi` not
                let (mut lo, mut hi) = if above {
                    (2.0 * nf, nf)
                } else {
                    (0.5 * nf, nf)
                };
                for _ in 0..200 {
                    let mid = 0.5 * (lo + hi);
                    if debye_valid(nf, mid) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                let x = hi;
                assert!(!debye_valid(nf, x) && x >= DEBYE_BELOW_MIN * nf);

                // ζ and t (see `olver`)
                let prec = 256;
                let n_r = rug::Float::with_val(prec, nf);
                let x_r = rug::Float::with_val(prec, x);
                let r2 = rug::Float::with_val(prec, n_r.square_ref()) - x_r.square();
                let y = r2.abs().sqrt() / &n_r;
                let c = if above {
                    rug::Float::with_val(prec, &y - y.clone().atan())
                } else {
                    rug::Float::with_val(prec, y.clone().atanh() - &y)
                };
                let zeta = (c * 3u32 / 2u32).pow(rug::Float::with_val(prec, 2) / 3u32);
                let zeta = if above { -zeta } else { zeta };
                let t = rug::Float::with_val(prec, n_r.pow(rug::Float::with_val(prec, 2) / 3u32))
                    * &zeta;
                let (t, zeta) = (t.to_f64(), zeta.to_f64());
                t_min = t_min.min(t);
                t_max = t_max.max(t);
                z_min = z_min.min(zeta);
                z_max = z_max.max(zeta);
            }
        }
        // The ranges documented in `olver`, within the ones of `AIRY_TBL` and
        // the polynomials of the coefficients
        assert!(-25.4 < t_min && t_max < 30.9);
        assert!(-0.218 < z_min && z_max < 0.304);
    }

    #[test]
    fn test_debye_below_min() {
        // From `ORDER_DEBYE_MIN`, the results with the arguments below
        // `DEBYE_BELOW_MIN * n` are zero (`jn`) or infinite (`yn`), and the
        // checks find it out (see `DEBYE_BELOW_MIN`). Their bounds get closer
        // to the limits as the argument increases (up to the order), and
        // farther from them, in proportion to the order, as it increases, so
        // the largest argument of each order is enough.
        let orders = (ORDER_DEBYE_MIN..=1 << 13).chain([1 << 20, (1 << 31) - 1, 1 << 31]);
        for n in orders {
            let x = f64::from_bits((DEBYE_BELOW_MIN * f64::from(n)).to_bits() - 1);
            assert!(jn_underflows(n, x, 1075), "n = {n}, x = {x}");
            assert!(yn_overflows(n, x, 1024), "n = {n}, x = {x}");
        }
    }

    #[test]
    fn test_olver() {
        // The results are compared before the final rounding, relative to
        // them up to the order and to the modulus beyond it, with a limit
        // somewhat above the documented error
        for n in [ORDER_DEBYE_MIN, 2000, 4096, 20_011] {
            let nf = f64::from(n);
            let below = debye_args(n, false).next().unwrap();
            let above = debye_args(n, true).next().unwrap();
            let spread = (0..=40).map(|i| below + (above - below) * (f64::from(i) / 40.0));
            for x in spread.chain([nf, f64::from_bits(above.to_bits() - 1)]) {
                if debye_valid(nf, x) {
                    continue;
                }
                let first = order_rug(n, x, false);
                let second = order_rug(n, x, true);
                let modulus = first.clone().hypot(&second);
                for (kind, exact) in [(false, first), (true, second)] {
                    let r = olver(n, x, kind);
                    let scale = if x <= nf {
                        exact.clone().abs()
                    } else {
                        modulus.clone()
                    };
                    let err = (float(r.hi()) + r.lo() - exact).abs() / scale;
                    let err = err.log2().to_f64();
                    assert!(
                        err < -65.0,
                        "n = {n}, x = {x:e}, second = {kind}, err = 2^{err}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_debye_turning() {
        // Close to `x = n`, where Olver's expansions are used, and around the
        // limits from which the Debye expansions are
        for n in [ORDER_DEBYE_MIN, 2000, 4096, 20_011] {
            let nf = f64::from(n);
            let below = debye_args(n, false).next().unwrap();
            let above = debye_args(n, true).next().unwrap();
            let xs = [
                below,
                f64::from_bits(below.to_bits() + 1),
                0.5 * (below + nf),
                nf - 1.0,
                nf,
                nf + 1.0,
                0.5 * (above + nf),
                f64::from_bits(above.to_bits() - 1),
                above,
            ];
            let spread = (1..40).map(|i| below + (above - below) * (f64::from(i) / 40.0));
            for x in xs.into_iter().chain(spread) {
                for second in [false, true] {
                    let actual = if second {
                        yn_debye(n, x)
                    } else {
                        jn_debye(n, x)
                    };
                    let exact = order_rug(n, x, second);
                    // In ULPs of the expected result. Up to the order, the
                    // error limit is the one of the other functions. Beyond
                    // it, the error can be larger by 2^-64 times the modulus.
                    let e = i32::from(exact.to_f64().exponent()).max(-1022) - 52;
                    let err = ((float(actual) - &exact).abs() >> e).to_f64();
                    let extra = if x <= nf {
                        0.0
                    } else {
                        let modulus = order_rug(n, x, false).hypot(&order_rug(n, x, true));
                        (modulus >> (64 + e)).to_f64()
                    };
                    let limit = 0.5 + f64::exp2i_fast(-10) + extra;
                    assert!(
                        err < limit,
                        "n = {n}, x = {x:e}, second = {second}, error = {err} ULP (limit = {limit})"
                    );
                }
            }
        }
    }

    /// Returns the terms of degree 1 to `DEBYE_TERMS` of `P(a, b)` (see
    /// `debye_poly`), with the coefficients that are used.
    fn debye_poly_rug(a: &rug::Float, b: &rug::Float) -> rug::Float {
        let mut acc = float(DEBYE_A1.hi()) + DEBYE_A1.lo();
        acc *= a;
        acc += rug::Float::with_val(PREC, b * DEBYE_B1.hi());
        let mut index = 0;
        for i in 0..=DEBYE_TERMS {
            for m in 0..=(DEBYE_TERMS - i) {
                let c = DEBYE_POLY[index];
                index += 1;
                let t = rug::Float::with_val(PREC, a.pow(m as u32))
                    * rug::Float::with_val(PREC, b.pow(i as u32));
                acc += t * c;
            }
        }
        assert_eq!(index, DEBYE_POLY.len());
        acc
    }

    #[test]
    fn test_debye_large_orders() {
        // The arithmetic of the expansions, compared with their evaluation
        // with `rug` (the functions of `rug` are too slow with these
        // orders). The phase of the functions of the first kind with these
        // orders is a large number, whose error is the largest with the
        // largest orders.
        for n in [(1 << 20) + 7, (1 << 31) - 1, 1 << 31] {
            let nf = f64::from(n);
            let n_r = float(nf);
            for x in debye_args(n, false) {
                let xr = float(x);
                let n2 = rug::Float::with_val(PREC, n_r.square_ref());
                let r = (n2 - rug::Float::with_val(PREC, xr.square_ref())).sqrt();
                let w = rug::Float::with_val(PREC, &r / &n_r);
                let b = rug::Float::with_val(PREC, r.recip_ref());
                let a = rug::Float::with_val(PREC, n_r.square_ref()) * b.clone().pow(3u32);
                let t = rug::Float::with_val(PREC, w.atanh_ref()) - &w;
                let core = rug::Float::with_val(PREC, &n_r * &t);
                for second in [false, true] {
                    if !debye_below_used(n, x, second) {
                        continue;
                    }
                    let e = debye_below(n, x, second);
                    let exact = if second {
                        let l = (rug::Float::with_val(PREC, &r * pi()) / 2u32).ln() / 2u32;
                        core.clone() - l + debye_poly_rug(&(-a.clone()), &(-b.clone()))
                    } else {
                        let l = (rug::Float::with_val(PREC, &r * pi()) * 2u32).ln() / 2u32;
                        -core.clone() - l + debye_poly_rug(&a, &b)
                    };
                    let err = (float(e.hi()) + e.lo() - exact).abs().log2().to_f64();
                    assert!(
                        err < -67.0,
                        "n = {n}, x = {x:e}, second = {second}, err = 2^{err}"
                    );
                }
            }
            for x in debye_args(n, true) {
                let prec = PREC + 1100;
                let pi = rug::Float::with_val(prec, rug::float::Constant::Pi);
                let xr = rug::Float::with_val(prec, x);
                let n_r = rug::Float::with_val(prec, nf);
                let x2 = rug::Float::with_val(prec, xr.square_ref());
                let r = (x2 - rug::Float::with_val(prec, n_r.square_ref())).sqrt();
                let v = rug::Float::with_val(prec, &r / &n_r);
                let b = rug::Float::with_val(prec, r.recip_ref());
                let a = rug::Float::with_val(prec, n_r.square_ref()) * b.clone().pow(3u32);
                let (p_re, p_im) = debye_poly_complex_rug(&a, &b);
                let theta = rug::Float::with_val(prec, &v - v.clone().atan()) * &n_r
                    - rug::Float::with_val(prec, &pi / 4u32)
                    + p_im;
                let modulus =
                    (rug::Float::with_val(prec, 2u32) / (r.clone() * &pi)).sqrt() * p_re.exp();

                let (m_hi, m_lo, k, c) = debye_above_parts(n, x);
                let err = rel_error_log2(&(float(m_hi) + m_lo), &modulus);
                assert!(err < -68.0, "n = {n}, x = {x:e}, modulus error = 2^{err}");
                let approx = rug::Float::with_val(prec, u32::from(k)) * &pi / 128u32 + c.bh + c.bl;
                // The difference is a multiple of `2 * π`
                let diff: rug::Float = (approx - theta) / (2u32 * pi.clone());
                let frac: rug::Float = diff.clone() - diff.round();
                let err = (frac * 2u32 * pi).abs().log2().to_f64();
                assert!(err < -66.0, "n = {n}, x = {x:e}, phase error = 2^{err}");
            }
        }
    }

    /// Returns `(re, im)` such that `re + i * im = P(i * a, -i * b)`, like
    /// `debye_poly_rug`.
    fn debye_poly_complex_rug(a: &rug::Float, b: &rug::Float) -> (rug::Float, rug::Float) {
        // The term a^m * b^i is multiplied by i^m * (-i)^i = i^(m + i) * (-1)^i
        let prec = a.prec();
        let mut re = rug::Float::new(prec);
        let mut im = (float(DEBYE_A1.hi()) + DEBYE_A1.lo()) * a;
        im -= rug::Float::with_val(prec, b * DEBYE_B1.hi());
        let mut index = 0;
        for i in 0..=DEBYE_TERMS {
            for m in 0..=(DEBYE_TERMS - i) {
                let c = DEBYE_POLY[index];
                index += 1;
                let t = rug::Float::with_val(prec, a.pow(m as u32))
                    * rug::Float::with_val(prec, b.pow(i as u32))
                    * c;
                // i^(m + i) * (-1)^i
                let k = m + i;
                let sign = if (k / 2 + i) % 2 == 0 { 1 } else { -1 };
                if k % 2 == 0 {
                    re += t * sign;
                } else {
                    im += t * sign;
                }
            }
        }
        (re, im)
    }

    /// Returns `a[n] / a[0]` and `a[n] / a[1]` for the backward recurrence of
    /// `jn_backward` that starts at the order `m`, without rounding errors.
    fn backward_ratios(n: u32, m: u32, x: f64) -> (rug::Float, rug::Float) {
        let h = float(2.0) / float(x);
        let (mut b, mut a) = (float(0.0), float(1.0));
        let mut an = float(0.0);
        for j in (1..=m).rev() {
            if j == n {
                an = a.clone();
            }
            let c = rug::Float::with_val(PREC, &h * j) * &a - &b;
            b = a;
            a = c;
        }
        (rug::Float::with_val(PREC, &an / &a), an / b)
    }

    /// Checks the truncation error of the backward recurrence of
    /// `jn_backward` when it starts at the order that `jn_backward_terms`
    /// gives, against a much later start.
    #[test]
    fn test_jn_backward_terms() {
        // Orders right below a power of 8 have the lowest threshold of
        // `q[k]` relative to their cube root
        for n in [2, 3, 7, 63, 511, 4095, 32_767, 262_143] {
            // The error is the largest right below the order
            let nf = f64::from(n);
            let below = f64::from_bits(nf.to_bits() - 1);
            for x in [below, 0.999 * nf, 0.9 * nf, 0.5 * nf, 0.01 * nf].map(f64::purify) {
                let k = jn_backward_terms(n, x);
                let (r0, r1) = backward_ratios(n, n + k, x);
                let (e0, e1) = backward_ratios(n, n + 3 * k + 50, x);
                // The one that is normalized with the largest of `j0(x)`
                // and `j1(x)`
                let err_log2 = if e0.clone().abs() <= e1.clone().abs() {
                    rel_error_log2(&r0, &e0)
                } else {
                    rel_error_log2(&r1, &e1)
                };
                assert!(err_log2 < -80.0, "n = {n}, x = {x:e}, error = 2^{err_log2}");
            }
        }
    }

    /// Returns arguments spread over `(0, n)`, from tiny ones to the largest
    /// one.
    fn below_order(n: u32) -> impl Iterator<Item = f64> {
        const RATIOS: [f64; 20] = [
            1e-300, 1e-100, 1e-30, 1e-6, 1e-3, 0.03, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
            0.95, 0.99, 0.999, 0.9999, 0.99999,
        ];
        let nf = f64::from(n);
        let below = f64::from_bits(nf.to_bits() - 1);
        RATIOS
            .into_iter()
            .map(move |z| (z * nf).purify())
            .chain([below])
    }

    /// Returns `log2(jn(n, x))`, for `n >= 1` and `0 < x < n`.
    fn jn_log2(n: u32, x: f64) -> f64 {
        let k = jn_backward_terms(n, x);
        let (r0, r1) = backward_ratios(n, n + 3 * k + 50, x);
        let (j0, j1) = (float(x).j0(), float(x).j1());
        // Normalized with the largest of `j0(x)` and `j1(x)`
        let j = if j0.clone().abs() >= j1.clone().abs() {
            j0 * r0
        } else {
            j1 * r1
        };
        j.log2().to_f64()
    }

    /// Returns `log2(-yn(n, x))`, for `n >= 2` and `0 < x < n`.
    fn yn_neg_log2(n: u32, x: f64) -> f64 {
        let xr = float(x);
        let (mut a, mut b) = (xr.clone().y0(), xr.clone().y1());
        for k in 1..n {
            let c = rug::Float::with_val(PREC, &b * (2 * k)) / &xr - &a;
            a = b;
            b = c;
        }
        (-b).log2().to_f64()
    }

    #[test]
    fn test_jn_log2_bound() {
        // Also with the largest orders, whose arguments close to the order
        // have the terms that cancel the most
        for n in [
            1,
            2,
            63,
            4095,
            262_143,
            10_000_000,
            i32::MAX as u32,
            1 << 31,
        ] {
            let nf = f64::from(n);
            for x in below_order(n).chain([nf]) {
                let z = float(x) / n;
                let w = (float(1.0) - z.clone().square()).sqrt();
                let ratio = z / (float(1.0) + &w);
                let exact = (ratio.log2() + w / float(2.0).ln()) * n;

                let err = (float(jn_log2_bound(n, x)) - &exact).abs().to_f64();
                let limit = nf * f64::exp2i_fast(-49) + exact.to_f64().abs() * f64::exp2i_fast(-50);
                assert!(err <= limit, "n = {n}, x = {x:e}, error = {err:e}");
            }
        }
    }

    #[test]
    fn test_jn_underflows() {
        for n in [1, 2, 3, 7, 63, 511, 4095, 32_767, 262_143] {
            // The bound exceeds the function by a factor of up to about
            // `sqrt(2 * π * n)`, and `jn_underflows` leaves a margin of one
            // more bit
            let nf = f64::from(n);
            let slack = 0.5 * crate::log2(2.0 * core::f64::consts::PI * nf) + 1.5;
            for x in below_order(n) {
                let exp = -jn_log2(n, x);
                // jn(n, x) >= 2^-e
                let e = crate::ceil(exp) as i32;
                assert!(!jn_underflows(n, x, e), "n = {n}, x = {x:e}, e = {e}");
                // jn(n, x) is far enough below 2^-e
                let e = crate::floor(exp - slack) as i32;
                assert!(jn_underflows(n, x, e), "n = {n}, x = {x:e}, e = {e}");
            }
        }
    }

    #[test]
    fn test_yn_overflows() {
        for n in [2, 3, 7, 63, 511, 4095, 32_767, 262_143] {
            // On top of the factor between the bound and the function of
            // `jn_underflows`, the term of the Wronskian that is left out is
            // the largest with the arguments close to the order, and both
            // factors together are up to about `n^(3/4)`
            let nf = f64::from(n);
            let slack = 0.75 * crate::log2(nf) + 4.5;
            for x in below_order(n).chain([nf - 1.0]) {
                let exp = yn_neg_log2(n, x);
                // yn(n, x) >= -2^e
                let e = crate::ceil(exp) as i32;
                assert!(!yn_overflows(n, x, e), "n = {n}, x = {x:e}, e = {e}");
                // yn(n, x) is far enough below -2^e. It is not checked
                // beyond `n - 1`, where the magnitude of the result is small.
                if x <= nf - 1.0 {
                    let e = crate::floor(exp - slack) as i32;
                    assert!(yn_overflows(n, x, e), "n = {n}, x = {x:e}, e = {e}");
                }
            }
        }
    }
}
