//! Gamma functions for `f64`.
//!
//! # `ln_gamma`
//!
//! For `x >= LN_GAMMA_STIRLING_MIN`, it uses the Stirling series (see
//! `stirling`):
//!
//! `ln(Γ(x)) = (x - 0.5) * ln(x) - x + ln(2 * π) / 2 + K1 / x + K3 / x^3 + ...`
//!
//! For smaller `x`, `x = n + f`, with an integer `n` and `|f| <= 1/2`, so
//!
//! `ln|Γ(x)| = ln(Γ(2 + f)) + ln(P)`, with `P = (2 + f) * ... * (x - 1)`,
//! when `n >= 3`, or
//!
//! `ln|Γ(x)| = ln(Γ(2 + f)) - ln|P|`, with `P = x * ... * (1 + f)`, when
//! `n <= 1`.
//!
//! `ln(Γ(2 + f))` is approximated with a polynomial, and `P` is calculated as
//! a sum of two `f64` (see `Prod`).
//!
//! For `x < LN_GAMMA_REFLECTION`, it uses the reflection formula
//!
//! `ln|Γ(x)| = ln(π) - ln|sin(π * x)| - ln(Γ(1 - x))`
//!
//! where `ln(Γ(1 - x))` is calculated with the Stirling series.
//!
//! The logarithms are calculated with `super::log::ln_accurate_parts` (except
//! in `ln_gamma_huge`, where the error of `super::log::ln_parts` is small
//! enough), so the error before the final rounding is less than about 2^-64.7,
//! plus about 2^-71 relative to the result (from the Stirling series). It is
//! small relative to the result except close to the zeros of `ln|Γ(x)|`.
//! Around 1 and 2, the relative errors of the terms are small. Around the
//! zeros with negative `x` (from -18 to -2), the terms nearly cancel, so the
//! results close to zero are calculated again with more accurate
//! approximations (see `ln_gamma_near_zero`), or with expansions around the
//! zeros for the smallest results (see `ln_gamma_near_root`).
//!
//! # `gamma`
//!
//! For `x >= GAMMA_STIRLING_MIN`, `Γ(x) = exp(ln(Γ(x)))`, where `ln(Γ(x))` is
//! calculated with an absolute error of less than about 2^-63.8 and the
//! exponential is evaluated like in `exp` (see `super::exp`).
//!
//! For `x <= GAMMA_REFLECTION`, it uses the reflection formula, also
//! evaluated as an exponential:
//!
//! `Γ(x) = ±exp(ln(π) - ln|sin(π * x)| - ln(Γ(1 - x)))`
//!
//! Otherwise, `x = n + f`, with an integer `n` and `|f| <= 1/2`, so
//!
//! `Γ(x) = Γ(3 + f) * P`, with `P = (3 + f) * ... * (x - 1)`, when `n >= 3`,
//! or
//!
//! `Γ(x) = Γ(3 + f) / P`, with `P = x * ... * (2 + f)`, when `n <= 2`.
//!
//! `Γ(3 + f)` is approximated with a polynomial, and `P` is calculated as a
//! sum of two `f64` (see `Prod`).
//!
//! The relative error before the final rounding is less than about 2^-63.

use super::exp::{Reduced as ExpReduced, exp_poly};
use super::f64x2::{F64x2, Split};
use super::log::{ln_accurate_parts, ln_accurate_sum_parts, ln_parts};
use super::round_i32;
use super::trigonometric::sinpi_f64x2;
use crate::traits::Float as _;

// GENERATE: consts F64x2 LN_PI
const LN_PI: F64x2 = F64x2::from_bits(0x3FF250D048E7A1BD, 0x3C67ABF2AD8D5088); // 1.144729885849400174143427351353e0

/// `ln_gamma` uses the Stirling series for `x >= LN_GAMMA_STIRLING_MIN`.
const LN_GAMMA_STIRLING_MIN: f64 = 4.0;

/// `gamma` uses the Stirling series for `x >= GAMMA_STIRLING_MIN`.
const GAMMA_STIRLING_MIN: f64 = 6.0;

/// `ln_gamma` uses the reflection formula for `x < LN_GAMMA_REFLECTION`.
const LN_GAMMA_REFLECTION: f64 = -4.0;

/// `gamma` uses the reflection formula for `x <= GAMMA_REFLECTION`.
const GAMMA_REFLECTION: f64 = -10.0;

// Coefficients of the Stirling series (see `stirling`), where `K0` and `K1`
// are exact (rounded to double-words), and the others are fitted for
// `x >= 4`, with an absolute error of 2^-68.5
// GENERATE: gamma_stirling_poly F64x2:2,f64 13 4
const K0: F64x2 = F64x2::from_bits(0x3FED67F1C864BEB5, 0xBC865B5A1B7FF5DF); // 9.189385332046727417803297364056e-1
const K1: F64x2 = F64x2::from_bits(0x3FB5555555555555, 0x3C55555555555555); // 8.333333333333333333333333333333e-2
const K3: f64 = f64::from_bits(0xBF66C16C16C16AB9); // -2.777777777777626e-3
const K5: f64 = f64::from_bits(0x3F4A01A019F1A8A5); // 7.936507935373439e-4
const K7: f64 = f64::from_bits(0xBF43813801D671FA); // -5.952380631018544e-4
const K9: f64 = f64::from_bits(0x3F4B9513E24C487C); // 8.417460527161548e-4
const K11: f64 = f64::from_bits(0xBF5F68E0A3C47A85); // -1.9170945880159506e-3
const K13: f64 = f64::from_bits(0x3F7A26DB1779B361); // 6.384712054210735e-3
const K15: f64 = f64::from_bits(0xBF9D32CB7FECD01F); // -2.8514079745221573e-2
const K17: f64 = f64::from_bits(0x3FC32A4B054BA70A); // 1.4972818144388383e-1
const K19: f64 = f64::from_bits(0xBFE8332759A1637D); // -7.56244349536345e-1
const K21: f64 = f64::from_bits(0x4006DF444DEA94FF); // 2.8590169989789724e0
const K23: f64 = f64::from_bits(0xC0162A6558BE58D6); // -5.541402231805668e0

impl crate::generic::Gamma for f64 {
    #[inline]
    fn gamma_finite(x: Self) -> Self {
        if x >= GAMMA_STIRLING_MIN {
            if x >= 172.0 {
                // Overflow
                f64::INFINITY
            } else {
                gamma_stirling(x)
            }
        } else if x > GAMMA_REFLECTION {
            gamma_small(x)
        } else if x > -185.0 {
            gamma_reflection(x)
        } else {
            // Underflow, with the sign of `Γ(x)`, which is negative between
            // `-2 * k - 1` and `-2 * k` (`x` is not an integer)
            let xi = (-x) as u64;
            0.0.set_sign((xi & 1) == 0)
        }
    }

    #[inline]
    fn ln_gamma_finite(x: Self) -> (Self, i8) {
        if x >= LN_GAMMA_STIRLING_MIN {
            (ln_gamma_stirling(x), 1)
        } else if x >= LN_GAMMA_REFLECTION {
            ln_gamma_small(x)
        } else {
            ln_gamma_reflection(x)
        }
    }
}

/// Returns `Γ(x)`, for `GAMMA_STIRLING_MIN <= x < 172`.
#[inline]
fn gamma_stirling(x: f64) -> f64 {
    // The absolute error of `ln(Γ(x))` is less than about 2^-63.8 (see
    // `stirling`: 2^-66, plus 2^-75 relative to (x - 0.5) * ln(x), which is
    // less than 883, and the error of the logarithm times `x - 0.5`, which is
    // as large), and so is the relative error of the result, plus the error of
    // the exponential (2^-64.5).
    let (l_hi, l_lo) = ln_accurate_parts(x, 0);
    let (hi, lo) = stirling(x, -0.5, l_hi, l_lo);
    let y = F64x2::fast_add11(hi, lo);
    let red = ExpReduced::exp_sum(y.hi(), y.lo());
    // It can overflow when x > 171.62
    red.eval_to_f64(exp_poly(red.r))
}

/// Returns `Γ(x)`, for `-185 < x <= GAMMA_REFLECTION`, where `x` is not an
/// integer.
#[inline(never)]
fn gamma_reflection(x: f64) -> f64 {
    // Γ(x) = ±exp(ln|Γ(x)|), where the absolute error of `ln|Γ(x)|` is less
    // than about 2^-63.4 (see `ln_gamma_reflection_parts`, with
    // (0.5 - x) * ln(-x) < 969), and so is the relative error of the result,
    // plus the error of the exponential (2^-64.5).
    let (hi, lo, negative) = ln_gamma_reflection_parts(x);
    // `hi` can be close to zero (when |Γ(x)| ~= 1), so the sum is exact
    let w = F64x2::add11(hi, lo);
    let red = ExpReduced::exp_sum(w.hi(), w.lo());
    // It can underflow when x < -170
    let r = red.eval_to_f64(exp_poly(red.r));
    if negative { -r } else { r }
}

/// Returns `Γ(x)`, for `GAMMA_REFLECTION < x < GAMMA_STIRLING_MIN`, where `x`
/// is not zero or a negative integer.
#[inline]
fn gamma_small(x: f64) -> f64 {
    if x.abs() < f64::exp2i_fast(-64) {
        // Γ(x) = 1 / x - γ + O(x), where `γ * |x| < 2^-64.8`, which is
        // negligible relative to `1 / x`
        return 1.0 / x;
    }

    // x = n + f, where `f` is exact
    let (nf, n) = round_i32(x);
    let f = x - nf;
    let (g_hi, g_lo) = gamma_poly(f);

    if n >= 3 {
        if n == 3 {
            return g_hi + g_lo;
        }
        // Γ(x) = Γ(3 + f) * (3 + f) * ... * (x - 1), where the factors are
        // exact
        let mut v = x - 1.0;
        let mut p = Prod::new(v, 0.0);
        for _ in 4..n {
            v -= 1.0;
            p = p.mul(v);
        }
        p.mul_to_f64(g_hi, g_lo)
    } else {
        let p = if n == 2 {
            // Γ(x) = Γ(3 + f) / (2 + f) = Γ(3 + f) / x
            Prod::new(x, 0.0)
        } else {
            // Γ(x) = Γ(3 + f) / ((1 + f) * (2 + f) * f * (f - 1) * ... * x),
            // where the factors after `(2 + f)` (when n <= 0) are exact
            let (e_hi, e_lo) = mul_1pf_2pf(f);
            let mut p = Prod::new(e_hi, e_lo);
            let mut v = f;
            for _ in n..1 {
                p = p.mul(v);
                v -= 1.0;
            }
            p
        };
        p.div_to_f64(g_hi, g_lo)
    }
}

/// Returns `ln(Γ(x))`, for `x >= LN_GAMMA_STIRLING_MIN`.
#[inline]
fn ln_gamma_stirling(x: f64) -> f64 {
    if x < f64::exp2i_fast(52) {
        // The absolute error is less than about 2^-66, plus about 2^-75
        // relative to the result (see `stirling`), which is greater than 1.79
        let (l_hi, l_lo) = ln_accurate_parts(x, 0);
        let (hi, lo) = stirling(x, -0.5, l_hi, l_lo);
        hi + lo
    } else {
        ln_gamma_huge(x)
    }
}

/// Returns `ln(Γ(x))`, for `x >= 2^52`, which overflows when
/// `x > ~2.5563 * 10^305`.
#[inline(never)]
fn ln_gamma_huge(x: f64) -> f64 {
    // ln(Γ(x)) = x * (ln(x) - 1) - 0.5 * ln(x) + K0 + O(1 / x), where the
    // last terms are less than 2^-46 relative to the first one, so their
    // rounding errors are negligible
    let (l_hi, l_lo) = ln_parts(x, 0);
    // exact (ln(x) > 36)
    let lm1 = l_hi - 1.0;
    // `lm1 = a + b` and `x = x1 + x2`, where `a` and `x1` have 26 bits, so
    // `x1 * a` and `x2 * a` are exact (or `x1 * a` overflows)
    let a = lm1.split_hi();
    let b = (lm1 - a) + l_lo;
    let x1 = x.split_hi();
    let x2 = x - x1;
    x1 * a + ((x2 * a + x * b) + (K0.hi() - 0.5 * l_hi))
}

/// Returns `(ln|Γ(x)|, sign(Γ(x)))`, for `x < LN_GAMMA_REFLECTION`, where `x`
/// is not an integer.
#[inline(never)]
fn ln_gamma_reflection(x: f64) -> (f64, i8) {
    let (hi, lo, negative) = ln_gamma_reflection_parts(x);
    let sign = if negative { -1 } else { 1 };
    let r = hi + lo;
    // The absolute error is less than about 2^-64.7 when `x > -20` (see
    // `ln_gamma_reflection_parts`, with (0.5 - x) * ln(-x) < 62), which is
    // not small enough relative to results close to zero (from -18 to -4).
    // With `x <= -20`, the magnitude of the result is at least 9.
    if r.abs() < 1.0 / 16.0 {
        (ln_gamma_near_zero(x), sign)
    } else {
        (r, sign)
    }
}

/// Returns `(hi, lo, negative)`, where `hi + lo ~= ln|Γ(x)|` (not normalized)
/// and `negative` tells whether `Γ(x)` is negative, for `-2^52 < x <= -4`,
/// where `x` is not an integer, with the reflection formula.
///
/// The absolute error is less than about 2^-65, plus about 2^-74 relative to
/// `(0.5 - x) * ln(-x)`:
/// * `ln(Γ(1 - x))`: about 2^-66, plus about 2^-75 relative to
///   `(0.5 - x) * ln(-x)`, and the error of the logarithm (2^-75) times
///   `0.5 - x` (see `stirling`).
/// * `ln|sin(π * x)|`: the relative error of the sine (about 2^-66), plus
///   2^-75 times the logarithm, whose magnitude is at most 34 (since
///   |sin(π * x)| >= sin(π * 2^-50)).
#[inline(always)]
fn ln_gamma_reflection_parts(x: f64) -> (f64, f64, bool) {
    // Γ(x) = π / (sin(π * x) * Γ(1 - x)), and with y = -x:
    // ln|Γ(x)| = ln(π) - ln|sin(π * x)| - ln(Γ(1 + y))
    let y = -x;
    let (l_hi, l_lo) = ln_accurate_parts(y, 0);
    let (r_hi, r_lo) = stirling(y, 0.5, l_hi, l_lo);
    let s = sinpi_f64x2(x);
    let negative = s.hi().is_sign_negative();
    let s = s.abs();
    let (ls_hi, ls_lo) = ln_accurate_sum_parts(s.hi(), s.lo());

    // The magnitudes can be close, and so can the result be to zero, so the
    // sums are exact
    let u = F64x2::add11(-r_hi, -ls_hi);
    let v = F64x2::add11(LN_PI.hi(), u.hi());
    let lo = ((v.lo() + u.lo()) + LN_PI.lo()) - (r_lo + ls_lo);
    (v.hi(), lo, negative)
}

/// Returns `(ln|Γ(x)|, sign(Γ(x)))`, for
/// `LN_GAMMA_REFLECTION <= x < LN_GAMMA_STIRLING_MIN`, where `x` is not zero
/// or a negative integer.
#[inline]
fn ln_gamma_small(x: f64) -> (f64, i8) {
    if x.abs() < f64::exp2i_fast(-60) {
        // ln|Γ(x)| = -ln|x| - γ * x + O(x^2), where `γ * |x| < 2^-60.8`,
        // which is negligible relative to `ln|x|` (greater than 41)
        let sign = if x.is_sign_negative() { -1 } else { 1 };
        let (y, edelta) = x.abs().normalize_arg();
        let (hi, lo) = ln_accurate_parts(y, edelta);
        return (-(hi + lo), sign);
    }

    // x = n + f, where `f` is exact
    let (nf, n) = round_i32(x);
    let f = x - nf;
    // |g| < 0.29
    let (g_hi, g_lo) = ln_gamma_poly(f);

    if n >= 2 {
        if n == 2 {
            return (g_hi + g_lo, 1);
        }
        // ln(Γ(x)) = ln(Γ(2 + f)) + ln((2 + f) * ... * (x - 1)), where the
        // factors are exact
        let mut v = x - 1.0;
        let (l_hi, l_lo) = if n == 3 {
            ln_accurate_parts(v, 0)
        } else {
            let mut p = Prod::new(v, 0.0);
            for _ in 3..n {
                v -= 1.0;
                p = p.mul(v);
            }
            let p = p.normalize();
            ln_accurate_sum_parts(p.hi(), p.lo())
        };
        // P >= 1.5 (when n = 3), so ln(P) >= 0.405 > |g|
        let s = F64x2::fast_add11(l_hi, g_hi);
        (s.hi() + (s.lo() + (l_lo + g_lo)), 1)
    } else if n == 1 {
        // ln(Γ(x)) = ln(Γ(2 + f)) - ln(1 + f) = ln(Γ(2 + f)) - ln(x)
        let (l_hi, l_lo) = ln_accurate_parts(x, 0);
        // |ln(1 + f)| > |ln(Γ(2 + f))| for |f| <= 1/2, and the result is
        // about -γ * f when `f` is close to zero, so the cancellation is
        // small
        let s = F64x2::fast_add11(-l_hi, g_hi);
        (s.hi() + (s.lo() + (g_lo - l_lo)), 1)
    } else {
        // ln|Γ(x)| = ln(Γ(2 + f)) - ln|f * (1 + f) * (f - 1) * ... * x|,
        // where the factors after `f * (1 + f)` (when n <= -1) are exact
        let (d_hi, d_lo) = mul_f_1pf(f);
        let mut p = Prod::new(d_hi, d_lo);
        let mut v = f;
        for _ in n..0 {
            v -= 1.0;
            p = p.mul(v);
        }
        let p = p.normalize();
        let sign = if p.hi().is_sign_negative() { -1 } else { 1 };
        let p = p.abs();
        let (l_hi, l_lo) = ln_accurate_sum_parts(p.hi(), p.lo());

        let r = sub_ln(g_hi, g_lo, l_hi, l_lo);
        // When the terms nearly cancel, the absolute error is less than about
        // 2^-67 (mostly from `ln(Γ(2 + f))`), which is not small enough
        // relative to results close to zero (from -4 to -2)
        if r.abs() < 1.0 / 128.0 {
            (ln_gamma_near_zero(x), sign)
        } else {
            (r, sign)
        }
    }
}

/// Returns `g_hi + g_lo - (l_hi + l_lo)`.
#[inline]
fn sub_ln(g_hi: f64, g_lo: f64, l_hi: f64, l_lo: f64) -> f64 {
    let s = F64x2::add11(g_hi, -l_hi);
    s.hi() + (s.lo() + (g_lo - l_lo))
}

/// Returns `ln|Γ(x)|`, for `-20 < x < -2`, when it is close to zero.
///
/// It is calculated like in `ln_gamma_small`, but with more accurate
/// approximations of `ln(Γ(2 + f))` and of the product, so the absolute error
/// is less than about 2^-76, which is small enough when the result is at
/// least 2^-16. Smaller results are calculated with an expansion around the
/// closest root (see `ln_gamma_near_root`).
#[cold]
#[inline(never)]
pub(crate) fn ln_gamma_near_zero(x: f64) -> f64 {
    // x = n + f, where `f` is exact
    let (nf, n) = round_i32(x);
    let f = x - nf;

    let (g_hi, g_lo) = ln_gamma_poly_accurate(f);

    // P = x * (x + 1) * ... * (1 + f), where the factors are exact (|x| > 2)
    let mut v = x;
    let mut p = F64x2::new1(x);
    for _ in n..1 {
        v += 1.0;
        p *= v;
    }
    let p = p.abs();
    let (l_hi, l_lo) = ln_accurate_sum_parts(p.hi(), p.lo());

    let r = sub_ln(g_hi, g_lo, l_hi, l_lo);
    if r.abs() >= f64::exp2i_fast(-16) {
        r
    } else {
        // `x` is always close enough to a root (see `LN_GAMMA_ROOTS`)
        ln_gamma_near_root(x).unwrap_or(r)
    }
}

/// Returns `ln|Γ(x)|` with an expansion around one of its negative roots (see
/// `LN_GAMMA_ROOTS`), if `x` is close enough to one of them.
///
/// Close to the roots, the result is the difference of two nearly equal
/// terms with the other methods, but the expansion has a relative error of
/// about 2^-64:
/// * Approximation: less than 2^-70.
/// * `x - x0` is calculated as a double-word (with `x0` as a sum of three
///   `f64`), and so is its product with the leading coefficient.
/// * The other terms are less than 2^-12.5 relative to the result.
#[cold]
fn ln_gamma_near_root(x: f64) -> Option<f64> {
    // The roots are sorted in decreasing order, so only the nearest ones
    // above and below `x` can be close enough
    let i = LN_GAMMA_ROOTS.partition_point(|root| root.x0[0] > x);
    for root in LN_GAMMA_ROOTS[i.saturating_sub(1)..].iter().take(2) {
        // When `x` is close enough to `x0`, both are in the same binade, so
        // the subtraction is exact
        let h = x - root.x0[0];
        if h.abs() <= root.radius {
            // x - x0 = d_hi + d_lo, where `h` is zero or a multiple of the
            // ULP of `x0[0]`, so it is larger than `x0[1]` in magnitude
            let d = F64x2::fast_add11(h, -root.x0[1]);
            let (d_hi, d_lo) = (d.hi(), d.lo() - root.x0[2]);

            // ln|Γ(x)| ~= d * k0 + d^2 * Q(d)
            let [k1, k2, k3, k4] = root.k;
            let q = k1 + d_hi * (k2 + d_hi * (k3 + d_hi * k4));
            let p = F64x2::mul11(d_hi, root.k0.hi());
            let lo = p.lo() + ((d_hi * root.k0.lo() + d_lo * root.k0.hi()) + (d_hi * d_hi) * q);
            return Some(p.hi() + lo);
        }
    }
    None
}

/// Returns `(hi, lo)` such that `hi + lo ~= ln(Γ(y + c + 1/2))`, where
/// `c = ±1/2` (so it is `ln(Γ(y))` or `ln(Γ(1 + y))`), for
/// `4 <= y < 2^52`, given `l_hi + l_lo ~= ln(y)` (with
/// `|l_lo| <= 2^-15 * |l_hi|`).
///
/// `ln(Γ(y + c + 1/2)) = (y + c) * ln(y) - y + K0 + K1 / y + K3 / y^3 + ...`
///
/// where the terms after `K1 / y` are approximated with a polynomial in
/// `1 / y^2`, with an absolute error of 2^-68.5.
///
/// The leading parts of the terms are added exactly, and the other parts are
/// less than 2^-14.5, or 2^-23 relative to `(y + c) * ln(y)`, so the absolute
/// error is less than about 2^-66, plus about 2^-75 relative to
/// `(y + c) * ln(y)` and the error of the logarithm times `y + c`.
///
/// `|lo| <= 2^-13 * |hi|`.
#[inline(always)]
fn stirling(y: f64, c: f64, l_hi: f64, l_lo: f64) -> (f64, f64) {
    // y + c = yc + yc_lo, where `yc_lo` is zero, except when the sum crosses
    // a power of two and its last bit is lost (when `c > 0`)
    let yc = (y + c).purify();
    let yc_lo = (y - yc) + c;
    // (y + c) * ln(y) = yc1 * a + (yc2 * a + yc * b + yc_lo * l_hi), where
    // `yc1` and `a` have 26 bits, so `yc1 * a` and `yc2 * a` are exact
    let yc1 = yc.split_hi();
    let yc2 = yc - yc1;
    let a = l_hi.split_hi();
    let b = (l_hi - a) + l_lo;

    // K3 / y^3 + K5 / y^5 + ...
    let t = 1.0 / y;
    let s = t * t;
    let s2 = s * s;
    let s4 = s2 * s2;
    let s8 = s4 * s4;
    let r = (((K3 + s * K5) + s2 * (K7 + s * K9)) + s4 * ((K11 + s * K13) + s2 * (K15 + s * K17)))
        + s8 * ((K19 + s * K21) + s2 * K23);
    let q3 = (t * s) * r;

    // K1 / y = q + (K1 - q * y) / y, where `q ~= K1 / y`
    //
    // q = q1 + q2 and y = y1 + y2, where `q1` and `y1` have 26 bits, so
    // `K1 - q * y` is calculated with the exact products `q1 * y1`, `q1 * y2`
    // and `q2 * y1` (`K1 - q1 * y1` is exact by Sterbenz lemma) and the
    // rounded `q2 * y2`, which is less than 2^-50 * K1
    let q = (K1.hi() * t).purify();
    let q1 = q.split_hi();
    let q2 = q - q1;
    let y1 = y.split_hi();
    let y2 = y - y1;
    let e = (((K1.hi() - q1 * y1) - q1 * y2) - q2 * y1) - q2 * y2;
    let q_lo = (e + K1.lo()) * t + q3;

    // yc1 * a >= 1.2 * y, and (y + c) * ln(y) - y is greater than K0 ~= 0.92
    // or in [0.5, 1) (when y < 4.2), so the sums with `fast_add11` are exact
    let s1 = F64x2::fast_add11(yc1 * a, -y);
    let s2 = F64x2::fast_add11(s1.hi(), K0.hi());
    let s3 = F64x2::fast_add11(s2.hi(), q);
    let lo =
        ((s1.lo() + s2.lo()) + s3.lo()) + (((yc2 * a + yc * b) + yc_lo * l_hi) + (K0.lo() + q_lo));
    (s3.hi(), lo)
}

/// Returns `(hi, lo)` such that `hi + lo ~= Γ(3 + f)`, for `|f| <= 1/2`, with
/// a relative error of about 2^-65.
///
/// `|lo| <= 2^-24 * |hi|`.
#[inline]
fn gamma_poly(f: f64) -> (f64, f64) {
    // GENERATE: gamma_poly F64x2:7,f64 17 3 -0.5 0.5
    const K0: F64x2 = F64x2::from_bits(0x4000000000000000, 0x3BDC830872F3B5F2); // 2.000000000000000000024150474365e0
    const K1: F64x2 = F64x2::from_bits(0x3FFD8773039049E7, 0x3C65669E453966C3); // 1.845568670196934278182285830791e0
    const K2: F64x2 = F64x2::from_bits(0x3FF3F18547938453, 0xBC647F2EC2EE9DEB); // 1.246464995951346512718587898809e0
    const K3: F64x2 = F64x2::from_bits(0x3FE2665A2BDCF38C, 0xBC887EFBCC2D3D2A); // 5.749941689206123493216013081455e-1
    const K4: F64x2 = F64x2::from_bits(0x3FCD73187D158329, 0xBC696C616D56D3E2); // 2.300749407541158476286906246152e-1
    const K5: F64x2 = F64x2::from_bits(0x3FB2DEFD42705C37, 0x3C532B928DC903C5); // 7.371504661601667867574312452513e-2
    const K6: F64x2 = F64x2::from_bits(0x3F9691F1CF9A570E, 0xBC1D73835B21E6A1); // 2.204110936744201708680002068834e-2
    const K7: f64 = f64::from_bits(0x3F76516EC8F77C26); // 5.448754076012674e-3
    const K8: f64 = f64::from_bits(0x3F56343551F0CC3C); // 1.3552208618079617e-3
    const K9: f64 = f64::from_bits(0x3F315A5DC0574561); // 2.647856602751766e-4
    const K10: f64 = f64::from_bits(0x3F100B444C2F8B5D); // 6.120304446670167e-5
    const K11: f64 = f64::from_bits(0x3EE1D66916B0A17C); // 8.505602371825604e-6
    const K12: f64 = f64::from_bits(0x3EC42F7EED7FBE10); // 2.406302751715784e-6
    const K13: f64 = f64::from_bits(0x3E77990463C8CE10); // 8.7908369424826e-8
    const K14: f64 = f64::from_bits(0x3E7E872F1540514D); // 1.1372589386909672e-7
    const K15: f64 = f64::from_bits(0xBE5128E8081B17E7); // -1.5981300289004878e-8
    const K16: f64 = f64::from_bits(0x3E44C413CF11CFB0); // 9.669888419704823e-9
    const K17: f64 = f64::from_bits(0xBE294A4E2172D65D); // -2.9441731375686135e-9

    // The terms of degree 7 and higher (less than 2^-14.5 relative to the
    // result) are evaluated with plain `f64` arithmetic
    let f2 = f * f;
    let f4 = f2 * f2;
    let t = ((K7 + f * K8) + f2 * (K9 + f * K10))
        + f4 * (((K11 + f * K12) + f2 * (K13 + f * K14)) + f4 * ((K15 + f * K16) + f2 * K17));

    let f = Split::new(f);
    let (h, l) = f.mul_add(t, 0.0, K6);
    let (h, l) = f.mul_add(h, l, K5);
    let (h, l) = f.mul_add(h, l, K4);
    let (h, l) = f.mul_add(h, l, K3);
    let (h, l) = f.mul_add(h, l, K2);
    let (h, l) = f.mul_add(h, l, K1);
    f.mul_add(h, l, K0)
}

/// Returns `(hi, lo)` such that `hi + lo ~= ln(Γ(2 + f))`, for `|f| <= 1/2`,
/// with a relative error of about 2^-64.5 and an absolute error of about
/// 2^-67.5.
///
/// `|lo| <= 2^-23 * |hi|`.
#[inline]
fn ln_gamma_poly(f: f64) -> (f64, f64) {
    // GENERATE: ln_gamma_poly F64x2:6,f64 22 2 -0.5 0.50001
    const K1: F64x2 = F64x2::from_bits(0x3FDB0EE6072093CE, 0x3C56CAA9EA01B471); // 4.227843350984671393927251703278e-1
    const K2: F64x2 = F64x2::from_bits(0x3FD4A34CC4A60FA6, 0x3C718E7C4CFD0C8E); // 3.224670334241132182607551977200e-1
    const K3: F64x2 = F64x2::from_bits(0xBFB13E001A557606, 0xBC5CD7AC92C1D9C6); // -6.735230105319809438634980441435e-2
    const K4: F64x2 = F64x2::from_bits(0x3F951322AC7D8481, 0x3C314251DF27F15E); // 2.058080842778454041278338187850e-2
    const K5: F64x2 = F64x2::from_bits(0xBF7E404FC218F67C, 0xBC1E406B26361B22); // -7.385551028674105782692753028334e-3
    const K6: F64x2 = F64x2::from_bits(0x3F67ADD6EADB7240, 0xBC08084C84AE205C); // 2.890510330742196269191173004102e-3
    const K7: f64 = f64::from_bits(0xBF538AC5C2BF05C5); // -1.192753911695697e-3
    const K8: f64 = f64::from_bits(0x3F40B36AF85FA4D1); // 5.096695247150067e-4
    const K9: f64 = f64::from_bits(0xBF2D3FD4C7F5E5D4); // -2.2315475869642712e-4
    const K10: f64 = f64::from_bits(0x3F1A127B11F22E4C); // 9.94575134304839e-5
    const K11: f64 = f64::from_bits(0xBF078DE59595698F); // -4.492623220171958e-5
    const K12: f64 = f64::from_bits(0x3EF580DC4D6EB43C); // 2.0507203625609103e-5
    const K13: f64 = f64::from_bits(0xBEE3CBD0A0BF0304); // -9.439540939017086e-6
    const K14: f64 = f64::from_bits(0x3ED25990E9DF977A); // 4.374949325458105e-6
    const K15: f64 = f64::from_bits(0xBEC11A57C530D9CD); // -2.038824769036462e-6
    const K16: f64 = f64::from_bits(0x3EB0043826B11799); // 9.546567081449322e-7
    const K17: f64 = f64::from_bits(0xBE9E45D3FBEA43E2); // -4.510993640095314e-7
    const K18: f64 = f64::from_bits(0x3E8CB5B1F5EE92A9); // 2.1390428181252712e-7
    const K19: f64 = f64::from_bits(0xBE7983BF72CB5300); // -9.504944093890693e-8
    const K20: f64 = f64::from_bits(0x3E675200A9DB543B); // 4.343748581456297e-8
    const K21: f64 = f64::from_bits(0xBE60CBED411DDB20); // -3.128608495790717e-8
    const K22: f64 = f64::from_bits(0x3E516D528A022C82); // 1.623019640953303e-8

    // The terms of degree 7 and higher (less than 2^-13 relative to the
    // result) are evaluated with plain `f64` arithmetic
    let f2 = f * f;
    let f4 = f2 * f2;
    let f8 = f4 * f4;
    let t = (((K7 + f * K8) + f2 * (K9 + f * K10)) + f4 * ((K11 + f * K12) + f2 * (K13 + f * K14)))
        + f8 * (((K15 + f * K16) + f2 * (K17 + f * K18))
            + f4 * ((K19 + f * K20) + f2 * (K21 + f * K22)));

    let f = Split::new(f);
    let (h, l) = f.mul_add(t, 0.0, K6);
    let (h, l) = f.mul_add(h, l, K5);
    let (h, l) = f.mul_add(h, l, K4);
    let (h, l) = f.mul_add(h, l, K3);
    let (h, l) = f.mul_add(h, l, K2);
    let (h, l) = f.mul_add(h, l, K1);
    f.mul(h, l)
}

/// A product `h + l` of factors of up to 53 bits, where `h` has at most 26
/// significant bits.
///
/// Each multiplication adds about 2^-24 to `|l / h|`, and a relative error of
/// about `|l / h| * 2^-52`.
#[derive(Copy, Clone)]
struct Prod {
    h: f64,
    l: f64,
}

impl Prod {
    /// Returns the product `hi + lo`, where `|lo| <= 2^-23 * |hi|`.
    #[inline]
    fn new(hi: f64, lo: f64) -> Self {
        let h = hi.split_hi();
        Self {
            h,
            l: (hi - h) + lo,
        }
    }

    /// Returns `self * v`.
    #[inline]
    fn mul(self, v: f64) -> Self {
        // v = v1 + v2, where `v1` has 26 bits, so `h * v1` and `h * v2` are
        // exact. `h * v1` is truncated to 26 bits (the new `h`), which keeps
        // the dependency chain of `h` short, and the remainder (less than
        // 2^-25 relative to it) is added to `l`.
        let v1 = v.split_hi();
        let v2 = v - v1;
        let p = self.h * v1;
        let h = p.split_hi();
        let l = ((p - h) + self.h * v2) + self.l * v;
        Self { h, l }
    }

    /// Returns the product normalized as a double-word.
    #[inline]
    fn normalize(self) -> F64x2 {
        F64x2::fast_add11(self.h, self.l)
    }

    /// Returns `(g_hi + g_lo) * self` rounded to `f64`, where
    /// `|g_lo| <= 2^-23 * |g_hi|`.
    #[inline]
    fn mul_to_f64(self, g_hi: f64, g_lo: f64) -> f64 {
        // g_hi = g1 + g2, where `g1` has 26 bits, so `g1 * h` and `g2 * h`
        // are exact
        let g1 = g_hi.split_hi();
        let g2 = g_hi - g1;
        g1 * self.h + (g2 * self.h + (g_hi * self.l + g_lo * (self.h + self.l)))
    }

    /// Returns `(g_hi + g_lo) / self` rounded to `f64`, where
    /// `|g_lo| <= 2^-23 * |g_hi|`.
    #[inline]
    fn div_to_f64(self, g_hi: f64, g_lo: f64) -> f64 {
        // q ~= g / p, and the remainder is calculated with `q = q1 + q2`,
        // where `q1` has 26 bits, so `q1 * h` and `q2 * h` are exact, and
        // `g_hi - q1 * h` is exact by Sterbenz lemma
        let d = self.h + self.l;
        let q = (g_hi / d).purify();
        let q1 = q.split_hi();
        let q2 = q - q1;
        let r = ((g_hi - q1 * self.h) - q2 * self.h) + (g_lo - q * self.l);
        q + r / d
    }
}

/// Returns `(hi, lo)` such that `hi + lo ~= f * (1 + f) = f + f^2`, for
/// `|f| <= 1/2`, with a relative error of about 2^-78.
///
/// `|lo| <= 2^-24 * |hi|`.
#[inline]
fn mul_f_1pf(f: f64) -> (f64, f64) {
    // f = f1 + f2, where `f1` has 26 bits, so `f1^2` is exact, and
    // f^2 = f1^2 + f2 * (f + f1)
    let f1 = f.split_hi();
    let f2 = f - f1;
    let s = F64x2::fast_add11(f, f1 * f1);
    (s.hi(), s.lo() + f2 * (f + f1))
}

/// Returns `(hi, lo)` such that
/// `hi + lo ~= (1 + f) * (2 + f) = 2 + 3 * f + f^2`, for `|f| <= 1/2`, with a
/// relative error of about 2^-78.
///
/// `|lo| <= 2^-24 * |hi|`.
#[inline]
fn mul_1pf_2pf(f: f64) -> (f64, f64) {
    // f = f1 + f2, where `f1` has 26 bits, so `3 * f1` (28 bits) and `f1^2`
    // are exact, and 3 * f + f^2 = 3 * f1 + f1^2 + f2 * (3 + f + f1)
    let f1 = f.split_hi();
    let f2 = f - f1;
    let s1 = F64x2::fast_add11(2.0, 3.0 * f1);
    // s1 >= 0.5 >= f1^2
    let s2 = F64x2::fast_add11(s1.hi(), f1 * f1);
    (s2.hi(), (s1.lo() + s2.lo()) + f2 * (3.0 + (f + f1)))
}

/// Returns `(hi, lo)` such that `hi + lo ~= ln(Γ(2 + f))`, for `|f| <= 1/2`,
/// like `ln_gamma_poly`, but with a relative error of about 2^-76.
///
/// `|lo| <= 2^-23 * |hi|`.
#[inline]
fn ln_gamma_poly_accurate(f: f64) -> (f64, f64) {
    // GENERATE: ln_gamma_poly F64x2:11,f64 26 2 -0.5 0.50001
    const K1: F64x2 = F64x2::from_bits(0x3FDB0EE6072093CE, 0x3C56CB90630CB51D); // 4.227843350984671393934877409352e-1
    const K2: F64x2 = F64x2::from_bits(0x3FD4A34CC4A60FA6, 0x3C71873E065D357D); // 3.224670334241132182362140608346e-1
    const K3: F64x2 = F64x2::from_bits(0xBFB13E001A557607, 0x3C5FB6D19F87FA62); // -6.735230105319809513301531490166e-2
    const K4: F64x2 = F64x2::from_bits(0x3F951322AC7D8483, 0x3C3AEF6C0C42ECD5); // 2.058080842778454787622409081228e-2
    const K5: F64x2 = F64x2::from_bits(0xBF7E404FC218F5F2, 0x3C1A7294E553F731); // -7.385551028673985318355199469977e-3
    const K6: F64x2 = F64x2::from_bits(0x3F67ADD6EADB6C31, 0xBC069D31D42CB6D9); // 2.890510330741523639774652619568e-3
    const K7: F64x2 = F64x2::from_bits(0xBF538AC5C2BF8DF3, 0x3BFCBFF575DC9BD1); // -1.192753911703256367804419556659e-3
    const K8: F64x2 = F64x2::from_bits(0x3F40B36AF8639627, 0xBBE5C45066B1EC2B); // 5.096695247430213797425690391245e-4
    const K9: F64x2 = F64x2::from_bits(0xBF2D3FD4C76D4E44, 0x3BC4D5A619AEAEB6); // -2.231547584537909048135672132060e-4
    const K10: F64x2 = F64x2::from_bits(0x3F1A127B0F18A18A, 0x3B96E8BDF4575A30); // 9.945751278251348044568505631641e-5
    const K11: F64x2 = F64x2::from_bits(0xBF078DE5BD6F867D, 0xBB938F51495D56A3); // -4.492623673236788447675615353489e-5
    const K12: f64 = f64::from_bits(0x3EF580DCEE24A019); // 2.0507212760950773e-5
    const K13: f64 = f64::from_bits(0xBEE3CBC967589C01); // -9.4394883759074e-6
    const K14: f64 = f64::from_bits(0x3ED2597A483BBD81); // 4.374866992881157e-6
    const K15: f64 = f64::from_bits(0xBEC11B2E122A9BCB); // -2.039214579713471e-6
    const K16: f64 = f64::from_bits(0x3EB0064AC56D0DAD); // 9.551393042741194e-7
    const K17: f64 = f64::from_bits(0xBE9E2629F11A61F8); // -4.492562633933671e-7
    const K18: f64 = f64::from_bits(0x3E8C7728DBB5DEFE); // 2.1208425758276032e-7
    const K19: f64 = f64::from_bits(0xBE7AF22907E70CD0); // -1.0038144358016404e-7
    const K20: f64 = f64::from_bits(0x3E6993F0170D8002); // 4.764251826968985e-8
    const K21: f64 = f64::from_bits(0xBE5895A76047D38E); // -2.289617919041097e-8
    const K22: f64 = f64::from_bits(0x3E47A2AE4EB379F2); // 1.1006124420206417e-8
    const K23: f64 = f64::from_bits(0xBE346BDA5E74A6BF); // -4.754704608584811e-9
    const K24: f64 = f64::from_bits(0x3E22949CC6E18B30); // 2.163056892497337e-9
    const K25: f64 = f64::from_bits(0xBE1E12ADE298DC51); // -1.7504769945952672e-9
    const K26: f64 = f64::from_bits(0x3E0F732C290BD5D7); // 9.153123416258258e-10

    let t = K12
        + horner!(
            f,
            f,
            [
                K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26
            ]
        );

    let f = Split::new(f);
    let (h, l) = f.mul_add(t, 0.0, K11);
    let (h, l) = f.mul_add(h, l, K10);
    let (h, l) = f.mul_add(h, l, K9);
    let (h, l) = f.mul_add(h, l, K8);
    let (h, l) = f.mul_add(h, l, K7);
    let (h, l) = f.mul_add(h, l, K6);
    let (h, l) = f.mul_add(h, l, K5);
    let (h, l) = f.mul_add(h, l, K4);
    let (h, l) = f.mul_add(h, l, K3);
    let (h, l) = f.mul_add(h, l, K2);
    let (h, l) = f.mul_add(h, l, K1);
    f.mul(h, l)
}

/// Expansion of `ln|Γ(x)|` around one of its negative roots `x0`.
struct LnGammaRoot {
    /// `x0 ~= x0[0] + x0[1] + x0[2]`
    x0: [f64; 3],
    /// The expansion is used when `|x - x0| <= radius`.
    radius: f64,
    /// `ln|Γ(x0 + h)| ~= h * (k0 + h * Q(h))`, with
    /// `Q(h) = k[0] + h * (k[1] + h * (k[2] + h * k[3]))`
    k0: F64x2,
    k: [f64; 4],
}

/// Expansions around the negative roots of `ln|Γ(x)|`, sorted in decreasing
/// order, from -2.457 to -13.
///
/// The radii are such that `|ln|Γ(x)|| >= 2^-15` at the limits of the ranges
/// (see the tests), so any argument where `ln_gamma_near_zero` finds a result
/// less than 2^-16 in magnitude is within one of them. The other roots (two
/// between each pair of consecutive integers until -18) are so close to the
/// poles that `|ln|Γ(x)|| > 2^-14` with the representable arguments around
/// them.
static LN_GAMMA_ROOTS: [LnGammaRoot; 23] = [
    ROOT_2_457,
    ROOT_2_747,
    ROOT_3_143,
    ROOT_3_955,
    ROOT_4_039,
    ROOT_4_991,
    ROOT_5_008,
    ROOT_5_998,
    ROOT_6_001,
    ROOT_6_999,
    ROOT_7_000,
    ROOT_7_999,
    ROOT_8_000,
    ROOT_8_999,
    ROOT_9_000,
    ROOT_9_999,
    ROOT_10_000,
    ROOT_10_999,
    ROOT_11_000,
    ROOT_11_999,
    ROOT_12_000,
    ROOT_12_999,
    ROOT_13_000,
];

// GENERATE: ln_gamma_root ROOT_2_457 -2.4570247382208006 -15 5
const ROOT_2_457: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC003A7FC9600F86C), // -2.4570247382208006e0
        f64::from_bits(0xBC855F64F98AF8D0), // -3.7075610815513266e-17
        f64::from_bits(0xB91C4B0CD201366A), // -1.3622663121726005e-33
    ],
    radius: f64::from_bits(0x3F00000000000000), // 3.0517578125e-5
    k0: F64x2::from_bits(0x3FF83FE966AF535F, 0xBC8775909A3D97DE), // 1.515603448021657321637058001907e0
    k: [
        f64::from_bits(0x40136EEBB002F55D), // 4.8583209516339965e0
        f64::from_bits(0x3FF694A6058A7858), // 1.41129114307798e0
        f64::from_bits(0x4021718D7D98DE35), // 8.721782612715382e0
        f64::from_bits(0x4017339FE06F9F7D), // 5.800414568722656e0
    ],
};

// GENERATE: ln_gamma_root ROOT_2_747 -2.7476826467274127 -15 5
const ROOT_2_747: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC005FB410A1BD901), // -2.7476826467274127e0
        f64::from_bits(0x3C9A19A96D2E6F85), // 9.055340329338315e-17
        f64::from_bits(0x393140B4FF4B7D60), // 3.322761057167369e-33
    ],
    radius: f64::from_bits(0x3F00000000000000), // 3.0517578125e-5
    k0: F64x2::from_bits(0xBFFEA12DA904B18C, 0xBC9220130F58898A), // -1.914350185611598816494731521973e0
    k: [
        f64::from_bits(0x4023267F3C265A52), // 9.575189475709667e0
        f64::from_bits(0xC034185AC30C8BF2), // -2.0095134916842603e1
        f64::from_bits(0x404F504AD31B3573), // 6.262728346661142e1
        f64::from_bits(0xC06858845C768BA6), // -1.9476615737107676e2
    ],
};

// GENERATE: ln_gamma_root ROOT_3_143 -3.1435808883499798 -17 5
const ROOT_3_143: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC009260DBC9E59AF), // -3.14358088834998e0
        f64::from_bits(0xBCAF717CD335A7B3), // -2.1818179852331714e-16
        f64::from_bits(0xB94D32A2A65BFD63), // -1.1246581285745781e-32
    ],
    radius: f64::from_bits(0x3EE0000000000000), // 7.62939453125e-6
    k0: F64x2::from_bits(0x401F20A65F2FAC55, 0xBCA1D258E4B05C0D), // 7.781884658131350872139552596618e0
    k: [
        f64::from_bits(0x4039D4D2977150EF), // 2.5831338372387957e1
        f64::from_bits(0x405C1137124D5C5B), // 1.12268986297176e2
        f64::from_bits(0x408267203E3130B5), // 5.888907436221112e2
        f64::from_bits(0x40A99A633920B431), // 3.2771937952251396e3
    ],
};

// GENERATE: ln_gamma_root ROOT_3_955 -3.9552942848585979 -19 5
const ROOT_3_955: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC00FA471547C2FE5), // -3.955294284858598e0
        f64::from_bits(0xBC770D4561291237), // -1.999428391746348e-17
        f64::from_bits(0x3909E6FADBBC171A), // 6.2357435447617e-34
    ],
    radius: f64::from_bits(0x3EC0000000000000), // 1.9073486328125e-6
    k0: F64x2::from_bits(0xC034B99D966C5647, 0x3CD9CBA2450AF735), // -2.072506084580370567830049982051e1
    k: [
        f64::from_bits(0x406F76DEAE0436BE), // 2.517146825868894e2
        f64::from_bits(0xC0AD25359D4B2F38), // -3.7306047156806126e3
        f64::from_bits(0x40EE8F829FDB342A), // 6.258808201370419e4
        f64::from_bits(0xC13116F780F5EEF1), // -1.1199915037526453e6
    ],
};

// GENERATE: ln_gamma_root ROOT_4_039 -4.0393618397405371 -19 5
const ROOT_4_039: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC010284E78599581), // -4.039361839740537e0
        f64::from_bits(0x3CAE78C1E9E43CFE), // 2.1143995503980602e-16
        f64::from_bits(0xB932AC17BFD6BE92), // -3.5961421111626576e-33
    ],
    radius: f64::from_bits(0x3EC0000000000000), // 1.9073486328125e-6
    k0: F64x2::from_bits(0x403ACA5CF4921642, 0x3CCA46A2E0D931BC), // 2.679048088614059326062263393759e1
    k: [
        f64::from_bits(0x40744415CD813F8E), // 3.2425532293784715e2
        f64::from_bits(0x40B559B11B2A9C7C), // 5.465691820777134e3
        f64::from_bits(0x40F96D18E2F09A44), // 1.0414555540523777e5
        f64::from_bits(0x4140261EB61E6D28), // 2.11666942280354e6
    ],
};

// GENERATE: ln_gamma_root ROOT_4_991 -4.9915446405600479 -21 5
const ROOT_4_991: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC013F7577A6EEAFD), // -4.991544640560048e0
        f64::from_bits(0x3CA5DE5EAB7F12CF), // 1.5174411760571722e-16
        f64::from_bits(0xB914075F5E0494A2), // -9.643515906617392e-34
    ],
    radius: f64::from_bits(0x3EA0000000000000), // 4.76837158203125e-7
    k0: F64x2::from_bits(0xC05D224A3EF9E41F, 0xBCF9BE272A13FF57), // -1.165357816162436307966887795883e2
    k: [
        f64::from_bits(0x40BB533C678A3956), // 6.995235954894064e3
        f64::from_bits(0xC120D3F7FEE65D34), // -5.514199978512884e5
        f64::from_bits(0x418752A6F6B5A2C1), // 4.8911582838689335e7
        f64::from_bits(0xC1F13D5D173062C7), // -4.627747187024116e9
    ],
};

// GENERATE: ln_gamma_root ROOT_5_008 -5.0082181683225935 -21 5
const ROOT_5_008: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC014086A57F0B6D9), // -5.0082181683225935e0
        f64::from_bits(0xBC895262B72CA9CA), // -4.3926353491015815e-17
        f64::from_bits(0xB92BD98D5E0861AA), // -2.68183947324466e-33
    ],
    radius: f64::from_bits(0x3EA0000000000000), // 4.76837158203125e-7
    k0: F64x2::from_bits(0x405ED72E0829AE02, 0xBCDFDC1859AD47E5), // 1.233621845633533933373067983332e2
    k: [
        f64::from_bits(0x40BCECC32EC22F9B), // 7.404762432228682e3
        f64::from_bits(0x412253D8563F7264), // 6.00556168452811e5
        f64::from_bits(0x418A225DF41555E8), // 5.480748651041776e7
        f64::from_bits(0x41F3E01774A25037), // 5.335250762144584e9
    ],
};

// GENERATE: ln_gamma_root ROOT_5_998 -5.9986074800808753 -24 5
const ROOT_5_998: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC017FE92F591F40D), // -5.998607480080875e0
        f64::from_bits(0xBCB7DD4ED62CBD32), // -3.311862478893795e-16
        f64::from_bits(0x3932071C071A2146), // 3.4720224807210337e-33
    ],
    radius: f64::from_bits(0x3E70000000000000), // 5.960464477539063e-8
    k0: F64x2::from_bits(0xC08661F6A43A5E12, 0xBD20C437B83BCB08), // -7.162454304275472920808369085228e2
    k: [
        f64::from_bits(0x410F79DCB794F26F), // 2.5785158963956262e5
        f64::from_bits(0xC19D6E8088A19FFE), // -1.2344528215783688e8
        f64::from_bits(0x422EF5D309A701E8), // 6.648663368350372e10
        f64::from_bits(0xC2C15EA6B1399FEB), // -3.8196442460991836e13
    ],
};

// GENERATE: ln_gamma_root ROOT_6_001 -6.0013852944531552 -24 5
const ROOT_6_001: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC018016B25897C8D), // -6.001385294453155e0
        f64::from_bits(0x3C927E0F49A4BA72), // 6.415847287933042e-17
        f64::from_bits(0xB9172E1AB15A4D03), // -1.116080967205309e-33
    ],
    radius: f64::from_bits(0x3E70000000000000), // 5.960464477539063e-8
    k0: F64x2::from_bits(0x40869DE49E3AF2AA, 0x3D0954B69094658A), // 7.237366299252801239787083874464e2
    k: [
        f64::from_bits(0x410FCE23484CFD10), // 2.6054841030309396e5
        f64::from_bits(0x419DE503A3C37C40), // 1.2538698494090366e8
        f64::from_bits(0x422F9C7B5326FFE9), // 6.7884657043499825e10
        f64::from_bits(0x42C1D3D507A7FD1B), // 3.920301964082621e13
    ],
};

// GENERATE: ln_gamma_root ROOT_6_999 -6.9998015078906377 -27 5
const ROOT_6_999: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC01BFFCBF76B86F0), // -6.999801507890638e0
        f64::from_bits(0x3C6853B29347B806), // 1.0550130037400023e-17
        f64::from_bits(0xB900FA018051DD41), // -4.08696427365735e-34
    ],
    radius: f64::from_bits(0x3E40000000000000), // 7.450580596923828e-9
    k0: F64x2::from_bits(0xC0B3ABF7A5CEA91B, 0xBD58257B8ABD091F), // -5.035967373768125414506135236209e3
    k: [
        f64::from_bits(0x4168349A2550422D), // 1.269064116604718e7
        f64::from_bits(0xC223D91DADC98428), // -4.262348976475812e10
        f64::from_bits(0x42E24F3D63CB8812), // 1.6105233333356856e14
        f64::from_bits(0xC3A20427DF8CC302), // -6.49103216841163e17
    ],
};

// GENERATE: ln_gamma_root ROOT_7_000 -7.0001983334073250 -27 5
const ROOT_7_000: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC01C0033FDEDFE1F), // -7.000198333407325e0
        f64::from_bits(0x3CB20BB7D2324678), // 2.504354173632409e-16
        f64::from_bits(0x395F5536678D69D3), // 2.413795840298293e-32
    ],
    radius: f64::from_bits(0x3E40000000000000), // 7.450580596923828e-9
    k0: F64x2::from_bits(0x40B3B407AA387BD1, 0x3D4DA1E57343BA18), // 5.044029941110829196534889597092e3
    k: [
        f64::from_bits(0x41683E85DAAFBAD6), // 1.2710958833951395e7
        f64::from_bits(0x4223E552B5E3C226), // 4.2725890801879196e10
        f64::from_bits(0x42E25E42A4BB56E9), // 1.6156843454328728e14
        f64::from_bits(0x43A216A356797A44), // 6.517043669854541e17
    ],
};

// GENERATE: ln_gamma_root ROOT_7_999 -7.9999751970958206 -30 5
const ROOT_7_999: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC01FFFF97F8159CF), // -7.999975197095821e0
        f64::from_bits(0xBC8E54F415A91586), // -5.261737128572354e-17
        f64::from_bits(0xB9253A5D106F9A3E), // -2.0441803623138533e-33
    ],
    radius: f64::from_bits(0x3E10000000000000), // 9.313225746154785e-10
    k0: F64x2::from_bits(0xC0E3AF76FE4C2FAB, 0xBD77CC92F0B99EC6), // -4.031571854218778888702656397163e4
    k: [
        f64::from_bits(0x41C838E76CAAF123), // 8.127648893354839e8
        f64::from_bits(0xC2B3DE68B3256526), // -2.184596023843715e13
        f64::from_bits(0x43A255C052AF928F), // 6.605867635342888e17
        f64::from_bits(0xC4920C2A8489E996), // -2.1306755336455246e22
    ],
};

// GENERATE: ln_gamma_root ROOT_8_000 -8.0000248002706815 -30 5
const ROOT_8_000: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC02000034028B3F9), // -8.000024800270682e0
        f64::from_bits(0xBCBF60CB3CEC1CED), // -4.354586297860107e-16
        f64::from_bits(0x395EA26620D6B1CA), // 2.3599860861934562e-32
    ],
    radius: f64::from_bits(0x3E10000000000000), // 9.313225746154785e-10
    k0: F64x2::from_bits(0x40E3B088FED67718, 0xBD8505613BA2961F), // 4.032428110812435386139726380541e4
    k: [
        f64::from_bits(0x41C83A3893550EDC), // 8.12937510664516e8
        f64::from_bits(0x42B3E0078DB8ADA4), // 2.185292033041364e13
        f64::from_bits(0x43A257BEC9A2D787), // 6.608673944431338e17
        f64::from_bits(0x44920E9EA0E747BC), // 2.1318070374708774e22
    ],
};

// GENERATE: ln_gamma_root ROOT_8_999 -8.9999972442509773 -33 5
const ROOT_8_999: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC021FFFFA3884BD0), // -8.999997244250977e0
        f64::from_bits(0xBCAFF90C9D2AE925), // -2.2185620509727132e-16
        f64::from_bits(0x39430C0EFEF78C04), // 7.336677520259467e-33
    ],
    radius: f64::from_bits(0x3DE0000000000000), // 1.1641532182693481e-10
    k0: F64x2::from_bits(0xC11625EDFC63DB2F, 0x3DADA7FC3ED681DA), // -3.628754964746711058723154973106e5
    k: [
        f64::from_bits(0x422EA8C150480A7A), // 6.584013008402046e10
        f64::from_bits(0xC34C4B30E4BC55C1), // -1.592794543191949e16
        f64::from_bits(0x446D5FE469975F11), // 4.334922735717809e21
        f64::from_bits(0xC59043D21C442A85), // -1.2584376013594556e27
    ],
};

// GENERATE: ln_gamma_root ROOT_9_000 -9.0000027557148226 -33 5
const ROOT_9_000: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC02200005C7768FB), // -9.000002755714823e0
        f64::from_bits(0xBC9B5B610FFB70D4), // -9.491348611623208e-17
        f64::from_bits(0xB93DEB7AD09EC5EA), // -5.762352109706189e-33
    ],
    radius: f64::from_bits(0x3DE0000000000000), // 1.1641532182693481e-10
    k0: F64x2::from_bits(0x4116261203919440, 0x3D97D5E8272CFF9C), // 3.628845034850277060821948962531e5
    k: [
        f64::from_bits(0x422EA8F32FB7F586), // 6.584176431597954e10
        f64::from_bits(0x434C4B75EE68E2BA), // 1.5928538462012788e16
        f64::from_bits(0x446D6043FADB9E47), // 4.3351379344785434e21
        f64::from_bits(0x45904414419F9786), // 1.2585156926857864e27
    ],
};

// GENERATE: ln_gamma_root ROOT_9_999 -9.9999997244266297 -36 5
const ROOT_9_999: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC023FFFFF6C0D7C0), // -9.99999972442663e0
        f64::from_bits(0x3CC197CEA8C42D7D), // 4.883037618642443e-16
        f64::from_bits(0x3967072C5A292198), // 3.548028340923709e-32
    ],
    radius: f64::from_bits(0x3DB0000000000000), // 1.4551915228366852e-11
    k0: F64x2::from_bits(0xC14BAF7DA5F3795D, 0xBDE16A79518CAD5E), // -3.628795296492739241627813856807e6
    k: [
        f64::from_bits(0x4297F3E8791FA0D2), // 6.584086185960205e12
        f64::from_bits(0xC3EBA18BEFCAAA63), // -1.5928210978304629e19
        f64::from_bits(0x4541EDE14818CC84), // 4.335019100486485e25
        f64::from_bits(0xC698D1A9AC8FF0C8), // -1.2584725700869054e32
    ],
};

// GENERATE: ln_gamma_root ROOT_10_000 -10.0000002755730133 -36 5
const ROOT_10_000: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC0240000093F2777), // -1.0000000275573013e1
        f64::from_bits(0xBCB927B45D95E154), // -3.4909708332642057e-16
        f64::from_bits(0xB950780C21B6E452), // -1.2687206116063323e-32
    ],
    radius: f64::from_bits(0x3DB0000000000000), // 1.4551915228366852e-11
    k0: F64x2::from_bits(0x414BAF825A0C63B2, 0xBDC20323F100B504), // 3.628804703503095512110545000721e6
    k: [
        f64::from_bits(0x4297F3EC8AE05F2E), // 6.584103254039795e12
        f64::from_bits(0x43EBA192FA62A5C8), // 1.5928272914951848e19
        f64::from_bits(0x4541EDE75FA72282), // 4.335041576057105e25
        f64::from_bits(0x4698D1B437228344), // 1.2584807260219779e32
    ],
};

// GENERATE: ln_gamma_root ROOT_10_999 -10.9999999749478903 -40 5
const ROOT_10_999: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC025FFFFFF28CDD4), // -1.099999997494789e1
        f64::from_bits(0x3CAC9924A65AA486), // 1.9843998306985407e-16
        f64::from_bits(0xB938D05A4E458063), // -4.778979059252407e-33
    ],
    radius: f64::from_bits(0x3D70000000000000), // 9.094947017729282e-13
    k0: F64x2::from_bits(0xC18308A7D8EADB7C, 0x3E0A95A609876C9F), // -3.991679511467644494886729002007e7
    k: [
        f64::from_bits(0x4306A4938065BFD2), // 7.966753636167623e14
        f64::from_bits(0xC491F51F646980C5), // -2.120048613893936e22
        f64::from_bits(0x462005993B63749B), // 6.346916356909948e29
        f64::from_bits(0xC7AE7EE7DD82E21C), // -2.0267886032135812e37
    ],
};

// GENERATE: ln_gamma_root ROOT_11_000 -11.0000000250521062 -40 5
const ROOT_11_000: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC026000000D7322A), // -1.1000000025052106e1
        f64::from_bits(0xBCC8AECB2D37FF52), // -6.850849812286175e-16
        f64::from_bits(0xB92C97D472001B98), // -2.753413969507158e-33
    ],
    radius: f64::from_bits(0x3D70000000000000), // 9.094947017729282e-13
    k0: F64x2::from_bits(0x418308A827152450, 0xBE21233B372BFA68), // 3.991680488532316485214718768261e7
    k: [
        f64::from_bits(0x4306A493DD62402E), // 7.966755586232378e14
        f64::from_bits(0x4491F51FD307A7CD), // 2.1200493922973838e22
        f64::from_bits(0x46200599BEFB1B89), // 6.346919464047429e29
        f64::from_bits(0x47AE7EE9169971AA), // 2.0267898434834354e37
    ],
};

// GENERATE: ln_gamma_root ROOT_11_999 -11.9999999979123242 -43 5
const ROOT_11_999: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC027FFFFFFEE1127), // -1.1999999997912324e1
        f64::from_bits(0xBC9CE1F7906B30F5), // -1.0020693920103036e-16
        f64::from_bits(0x390B43A13E31B9DF), // 6.563612372549864e-34
    ],
    radius: f64::from_bits(0x3D40000000000000), // 1.1368683772161603e-13
    k0: F64x2::from_bits(0xC1BC8CFBFAF2B0C8, 0x3E47E94018C659D3), // -4.790015949480099566694655040309e8
    k: [
        f64::from_bits(0x43797926203E98AC), // 1.1472126519132435e17
        f64::from_bits(0xC53E4DA54EBC6DAC), // -3.6634446193922216e25
        f64::from_bits(0x470447163BBA4E0E), // 1.316096871962194e34
        f64::from_bits(0xC8CCF2769FE2A352), // -5.0433000359192426e42
    ],
};

// GENERATE: ln_gamma_root ROOT_12_000 -12.0000000020876758 -43 5
const ROOT_12_000: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC02800000011EED9), // -1.2000000002087676e1
        f64::from_bits(0x3CA19D5307E1FB5E), // 1.2222548112048185e-16
        f64::from_bits(0x3928F0DBE4153150), // 2.4017170001173477e-33
    ],
    radius: f64::from_bits(0x3D40000000000000), // 1.1368683772161603e-13
    k0: F64x2::from_bits(0x41BC8CFC050D4F38, 0xBE57EF0726D4B0CF), // 4.790016050519900099060471392974e8
    k: [
        f64::from_bits(0x4379792629426754), // 1.1472126761123565e17
        f64::from_bits(0x453E4DA55ED2869F), // 3.66344473530636e25
        f64::from_bits(0x470447164A14833E), // 1.316096927485252e34
        f64::from_bits(0x48CCF276B97EF7C6), // 5.0433003018755794e42
    ],
};

// GENERATE: ln_gamma_root ROOT_12_999 -12.9999999998394102 -47 5
const ROOT_12_999: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC029FFFFFFFE9EDC), // -1.299999999983941e1
        f64::from_bits(0x3CC84F40342D001C), // 6.747262033096337e-16
        f64::from_bits(0x39650556E5AEDE66), // 3.2387758664429733e-32
    ],
    radius: f64::from_bits(0x3D00000000000000), // 7.105427357601002e-15
    k0: F64x2::from_bits(0xC1F7328CBFACB4E5, 0xBE7EAE875D916AC1), // -6.227020794794163818216893872667e9
    k: [
        f64::from_bits(0x43F0D0FA2E06B2F1), // 1.9387894005607895e19
        f64::from_bits(0xC5F04105BEC453B2), // -8.048587946043371e28
        f64::from_bits(0x47F1AC9DD47DC9F8), // 3.7589043458591325e38
        f64::from_bits(0xC9F47FFB07535211), // -1.8725420437310025e48
    ],
};

// GENERATE: ln_gamma_root ROOT_13_000 -13.0000000001605898 -47 5
const ROOT_13_000: LnGammaRoot = LnGammaRoot {
    x0: [
        f64::from_bits(0xC02A000000016124), // -1.300000000016059e1
        f64::from_bits(0xBCC84E03341EE8DD), // -6.745919484964342e-16
        f64::from_bits(0x396F8391FEF50BD4), // 4.8554922539526397e-32
    ],
    radius: f64::from_bits(0x3D00000000000000), // 7.105427357601002e-15
    k0: F64x2::from_bits(0x41F7328CC0534B1B, 0xBE7F63C3A52BDF72), // 6.227020805205836179145782966731e9
    k: [
        f64::from_bits(0x43F0D0FA2E7F760F), // 1.9387894038024745e19
        f64::from_bits(0x45F04105BF7369B6), // 8.048587966229412e28
        f64::from_bits(0x47F1AC9DD57BA2C6), // 3.758904358429022e38
        f64::from_bits(0x49F47FFB08C35BB7), // 1.8725420515582985e48
    ],
};

#[cfg(test)]
mod tests {
    use super::LN_GAMMA_ROOTS;

    #[test]
    fn test_ln_gamma_roots() {
        // `ln_gamma_near_root` requires the roots to be sorted in decreasing
        // order, with non-overlapping ranges
        for pair in LN_GAMMA_ROOTS.windows(2) {
            assert!(pair[0].x0[0] - pair[0].radius > pair[1].x0[0] + pair[1].radius);
        }

        // `ln_gamma_near_zero` requires `|ln|Γ(x)|| >= 2^-15` at the limits of
        // the ranges (greater than 2^-16 plus its error)
        for root in LN_GAMMA_ROOTS.iter() {
            for x in [root.x0[0] - root.radius, root.x0[0] + root.radius] {
                let y = rug::Float::with_val(256, x).ln_abs_gamma().0;
                assert!(y.abs() >= rug::Float::with_val(256, rug::Float::i_exp(1, -15)));
            }
        }
    }

    #[test]
    fn test_ln_gamma_roots_complete() {
        // `ln_gamma_near_zero` also requires an expansion around every zero
        // of `ln|Γ(x)|` with representable arguments where `|ln|Γ(x)||` is
        // less than 2^-16 (plus its error). `|ln|Γ(x)||` increases away from
        // the zeros, so, with the requirement checked in
        // `test_ln_gamma_roots`, it is enough to check the arguments closest
        // to each zero.
        const PREC: u32 = 256;
        let ln_gamma = |x: &rug::Float| rug::Float::with_val(PREC, x).ln_abs_gamma().0;
        let limit = rug::Float::with_val(PREC, rug::Float::i_exp(1, -15));
        let bisect = |mut pos: rug::Float, mut neg: rug::Float| {
            for _ in 0..200 {
                let m = rug::Float::with_val(PREC, &pos + &neg) / 2;
                if ln_gamma(&m).is_sign_negative() {
                    neg = m;
                } else {
                    pos = m;
                }
            }
            pos
        };

        // Two zeros between each pair of consecutive integers from -20 to -2
        // (where `ln_gamma_near_zero` can be used), on both sides of the
        // minimum of `|Γ(x)|` (where ψ(x) = 0)
        for n in 2..20 {
            let mut a = rug::Float::with_val(PREC, -n - 1);
            let mut b = rug::Float::with_val(PREC, -n);
            for _ in 0..200 {
                let m = rug::Float::with_val(PREC, &a + &b) / 2;
                if rug::Float::with_val(PREC, &m).digamma().is_sign_negative() {
                    a = m;
                } else {
                    b = m;
                }
            }
            let left = bisect(rug::Float::with_val(PREC, -n - 1), a.clone());
            let right = bisect(rug::Float::with_val(PREC, -n), a);

            for zero in [left, right] {
                let x0 = zero.to_f64();
                for k in -8..=8 {
                    let x = f64::from_bits(x0.to_bits().wrapping_add_signed(k));
                    if x.fract() == 0.0 || ln_gamma(&rug::Float::with_val(PREC, x)).abs() >= limit {
                        continue;
                    }
                    assert!(
                        LN_GAMMA_ROOTS
                            .iter()
                            .any(|root| (x - root.x0[0]).abs() <= root.radius),
                        "no expansion around {x:e}",
                    );
                }
            }
        }
    }
}
