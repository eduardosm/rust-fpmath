mod cbrt;
mod exp;
mod frexp;
mod gamma;
mod hyperbolic;
mod hypot;
mod inv_hyperbolic;
mod inv_trigonometric;
mod log;
mod pow;
mod round;
mod scalbn;
mod sqrt;
mod trigonometric;
mod values;

const MIN_MAX_ERROR: f64 = 0.4999;
const ERROR_LIMIT: f64 = 0.51;

const RUG_PREC: u32 = 53 + 20;

fn check_result(
    input: impl std::fmt::Debug,
    actual: f64,
    expected: rug::Float,
    max_error: &mut f64,
) {
    let expected_is_neg = expected.is_sign_negative();
    let err = calc_error_ulp(actual, expected);
    *max_error = max_error.max(err);

    assert!(err < ERROR_LIMIT, "input = {input:?}, error = {err} ULP");
    if !actual.is_nan() {
        assert_eq!(
            actual.is_sign_negative(),
            expected_is_neg,
            "input = {input:?}, sign mismatch",
        );
    }
}

/// Calculates the error of `actual` in ULPs with respect to the reference
/// result `expected`.
///
/// The distance is measured with the local spacing of the values: each
/// interval between two consecutive values counts as one ULP, and the
/// position of `expected` inside its interval gives the fractional part. This
/// is `|actual - expected| / ulp(expected)` when both are in the same binade,
/// but it also works across powers of two, where the spacing changes: the
/// value before `2^k` is one ULP away from `2^k`, not half an ULP. So an error
/// smaller than 0.5 ULP means that the result is correctly rounded, and an
/// error of `0.5 + d` ULP means that it is rounded to the wrong side of a
/// midpoint that is `d` ULP away from the exact result.
///
/// Results that overflow are rounded to infinity, which is treated as the
/// value after `f64::MAX` (2^1024). So exact results in `(MAX, MAX + ulp/2)`
/// must be rounded to `MAX`.
fn calc_error_ulp(actual: f64, expected: rug::Float) -> f64 {
    let actual = purify(actual);
    let overflow: rug::Float = rug::Float::with_val(53, 1) << 1024;

    if expected.is_nan() {
        if actual.is_nan() { 0.0 } else { f64::INFINITY }
    } else if actual.is_nan() {
        f64::INFINITY
    } else if expected.clone().abs() >= overflow {
        if actual.is_infinite() && actual.is_sign_negative() == expected.is_sign_negative() {
            0.0
        } else {
            f64::INFINITY
        }
    } else {
        // Value with ordinal `ord`, with infinity at 2^1024.
        let value = |ord: i64| {
            let x = from_ordinal(ord);
            if x.is_infinite() {
                rug::Float::with_val(53, &overflow * x.signum())
            } else {
                rug::Float::with_val(53, x)
            }
        };
        // `expected` is between the consecutive values with ordinals `lo` and
        // `lo + 1`.
        let lo = ordinal(expected.to_f64_round(rug::float::Round::Down));
        let (lo_value, hi_value) = (value(lo), value(lo + 1));
        let prec = expected.prec() + 64;
        let frac = rug::Float::with_val(prec, &expected - &lo_value)
            / rug::Float::with_val(prec, &hi_value - &lo_value);
        ((i128::from(ordinal(actual)) - i128::from(lo)) as f64 - frac.to_f64()).abs()
    }
}

/// Position of `x` in the ordered sequence of values (`+0` and `-0` are
/// both mapped to zero).
fn ordinal(x: f64) -> i64 {
    let bits = x.to_bits();
    let mag = (bits & !0x8000_0000_0000_0000) as i64;
    if (bits & 0x8000_0000_0000_0000) != 0 {
        -mag
    } else {
        mag
    }
}

/// Inverse of [`ordinal`] (zero is mapped to `+0`).
fn from_ordinal(ord: i64) -> f64 {
    if ord < 0 {
        f64::from_bits(ord.unsigned_abs() | 0x8000_0000_0000_0000)
    } else {
        f64::from_bits(ord as u64)
    }
}

fn exponent(x: f64) -> i32 {
    let biased = ((x.to_bits() & !0x8000_0000_0000_0000) >> 52) as i32;
    biased - 1023
}

// Workaround X87 compiler bugs
fn purify(x: f64) -> f64 {
    std::hint::black_box(x)
}

impl crate::ResultEq for f64 {
    fn result_eq(&self, other: &Self) -> bool {
        if self.is_nan() && other.is_nan() {
            // NaNs with different payloads or signs are considered equal.
            true
        } else {
            self.to_bits() == other.to_bits()
        }
    }
}

#[test]
fn test_calc_error_ulp() {
    let exact = |x: f64| rug::Float::with_val(53, x);
    let near = |x: f64, ulps: f64| {
        // `x + ulps * ulp(x)`, with the ULP of the binade of `x`
        let ulp = rug::Float::with_val(53, 1) << (exponent(x) - 52);
        rug::Float::with_val(128, x) + ulp * ulps
    };
    let pred = |x: f64| from_ordinal(ordinal(x) - 1);
    let succ = |x: f64| from_ordinal(ordinal(x) + 1);

    // Same binade
    assert_eq!(calc_error_ulp(1.5, near(1.5, 0.375)), 0.375);
    assert_eq!(calc_error_ulp(succ(1.5), near(1.5, 0.375)), 0.625);
    // The values around a power of two are one ULP away from it, although
    // the spacing is halved below it.
    assert_eq!(calc_error_ulp(pred(1.0), exact(1.0)), 1.0);
    assert_eq!(calc_error_ulp(succ(1.0), exact(1.0)), 1.0);
    assert_eq!(calc_error_ulp(-pred(1.0), exact(-1.0)), 1.0);
    // Just above and just below a power of two
    assert_eq!(calc_error_ulp(1.0, near(1.0, 0.25)), 0.25);
    assert_eq!(calc_error_ulp(pred(1.0), near(1.0, 0.25)), 1.25);
    assert_eq!(calc_error_ulp(1.0, near(pred(1.0), 0.75)), 0.25);
    assert_eq!(calc_error_ulp(pred(1.0), near(pred(1.0), 0.75)), 0.75);
    // Infinity is the value after `MAX`
    assert_eq!(calc_error_ulp(f64::MAX, near(f64::MAX, 0.25)), 0.25);
    assert_eq!(calc_error_ulp(f64::INFINITY, near(f64::MAX, 0.25)), 0.75);
    assert_eq!(calc_error_ulp(-f64::INFINITY, near(-f64::MAX, -0.75)), 0.25);
    assert_eq!(calc_error_ulp(f64::INFINITY, exact(f64::INFINITY)), 0.0);
    // Subnormals and zero
    let min = f64::from_bits(1);
    assert_eq!(calc_error_ulp(min, exact(0.0)), 1.0);
    assert_eq!(calc_error_ulp(-min, exact(min)), 2.0);
    assert_eq!(calc_error_ulp(0.0, rug::Float::with_val(53, min) / 4), 0.25);
    assert_eq!(
        calc_error_ulp(f64::MIN_POSITIVE, exact(pred(f64::MIN_POSITIVE))),
        1.0
    );

    // Exact ties, zeros and NaN
    assert_eq!(calc_error_ulp(1.5, near(1.5, 0.5)), 0.5);
    assert_eq!(calc_error_ulp(succ(1.5), near(1.5, 0.5)), 0.5);
    assert_eq!(calc_error_ulp(0.0, exact(-0.0)), 0.0);
    assert_eq!(calc_error_ulp(-0.0, exact(0.0)), 0.0);
    assert_eq!(calc_error_ulp(f64::NAN, exact(f64::NAN)), 0.0);
    assert_eq!(calc_error_ulp(1.0, exact(f64::NAN)), f64::INFINITY);
    assert_eq!(calc_error_ulp(f64::NAN, exact(1.0)), f64::INFINITY);

    // In `[MAX + ulp/2, 2^1024)`, negative next to a power of two, and tiny
    // and negative (see the `f32` version, which rounds `expected` to
    // nearest instead of down).
    assert_eq!(calc_error_ulp(f64::MAX, near(f64::MAX, 0.75)), 0.75);
    assert_eq!(calc_error_ulp(f64::INFINITY, near(f64::MAX, 0.75)), 0.25);
    assert_eq!(calc_error_ulp(-1.0, near(-1.0, -0.25)), 0.25);
    assert_eq!(calc_error_ulp(-succ(1.0), near(-1.0, -0.25)), 0.75);
    assert_eq!(calc_error_ulp(-f64::MAX, near(-f64::MAX, -0.25)), 0.25);
    let tiny_neg: rug::Float = rug::Float::with_val(53, -min) / 4;
    assert_eq!(calc_error_ulp(-0.0, tiny_neg.clone()), 0.25);
    assert_eq!(calc_error_ulp(-min, tiny_neg), 0.75);
}
