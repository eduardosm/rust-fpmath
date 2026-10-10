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

const MIN_MAX_ERROR: f32 = 0.4999;
const ERROR_LIMIT: f32 = 0.5 + 1.0 / 2048.0;

fn check_result(input: impl std::fmt::Debug, actual: f32, expected: f64, max_error: &mut f32) {
    let err = calc_error_ulp(actual, expected);
    *max_error = max_error.max(err);

    assert!(err < ERROR_LIMIT, "input = {input:?}, error = {err} ULP");
    if !actual.is_nan() {
        assert_eq!(
            actual.is_sign_negative(),
            expected.is_sign_negative(),
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
/// value after `f32::MAX` (2^128). So exact results in `(MAX, MAX + ulp/2)`
/// must be rounded to `MAX`.
fn calc_error_ulp(actual: f32, expected: f64) -> f32 {
    let actual = purify(actual);
    let overflow = f64::from_bits((1023 + 128) << 52);

    if expected.is_nan() {
        if actual.is_nan() { 0.0 } else { f32::INFINITY }
    } else if actual.is_nan() {
        f32::INFINITY
    } else if expected.abs() >= overflow {
        if actual.is_infinite() && actual.is_sign_negative() == expected.is_sign_negative() {
            0.0
        } else {
            f32::INFINITY
        }
    } else {
        // Value with ordinal `ord`, with infinity at 2^128.
        let value = |ord: i32| {
            let x = f64::from(from_ordinal(ord));
            if x.is_infinite() {
                overflow.copysign(x)
            } else {
                x
            }
        };
        // `expected` is between the consecutive values with ordinals `lo` and
        // `lo + 1`.
        let mut lo = ordinal(purify(expected as f32));
        if value(lo) > expected {
            lo -= 1;
        }
        let (lo_value, hi_value) = (value(lo), value(lo + 1));
        // Exact: `lo_value` and `expected` are close, and the divisor is a
        // power of two.
        let frac = (expected - lo_value) / (hi_value - lo_value);
        ((i64::from(ordinal(actual)) - i64::from(lo)) as f64 - frac).abs() as f32
    }
}

/// Position of `x` in the ordered sequence of values (`+0` and `-0` are
/// both mapped to zero).
fn ordinal(x: f32) -> i32 {
    let bits = x.to_bits();
    let mag = (bits & !0x8000_0000) as i32;
    if (bits & 0x8000_0000) != 0 { -mag } else { mag }
}

/// Inverse of [`ordinal`] (zero is mapped to `+0`).
fn from_ordinal(ord: i32) -> f32 {
    if ord < 0 {
        f32::from_bits(ord.unsigned_abs() | 0x8000_0000)
    } else {
        f32::from_bits(ord as u32)
    }
}

fn exponent(x: f32) -> i32 {
    let biased = ((x.to_bits() & !0x8000_0000) >> 23) as i32;
    biased - 127
}

// Workaround X87 compiler bugs
fn purify(x: f32) -> f32 {
    std::hint::black_box(x)
}

impl crate::ResultEq for f32 {
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
    let near = |x: f32, ulps: f64| {
        // `x + ulps * ulp(x)`, with the ULP of the binade of `x`
        let ulp = f64::from_bits(((1023 + exponent(x) - 23) as u64) << 52);
        f64::from(x) + ulps * ulp
    };
    let pred = |x: f32| from_ordinal(ordinal(x) - 1);
    let succ = |x: f32| from_ordinal(ordinal(x) + 1);

    // Same binade
    assert_eq!(calc_error_ulp(1.5, near(1.5, 0.375)), 0.375);
    assert_eq!(calc_error_ulp(succ(1.5), near(1.5, 0.375)), 0.625);
    // The values around a power of two are one ULP away from it, although
    // the spacing is halved below it.
    assert_eq!(calc_error_ulp(pred(1.0), 1.0), 1.0);
    assert_eq!(calc_error_ulp(succ(1.0), 1.0), 1.0);
    assert_eq!(calc_error_ulp(-pred(1.0), -1.0), 1.0);
    // Just above and just below a power of two
    assert_eq!(calc_error_ulp(1.0, near(1.0, 0.25)), 0.25);
    assert_eq!(calc_error_ulp(pred(1.0), near(1.0, 0.25)), 1.25);
    assert_eq!(calc_error_ulp(1.0, near(pred(1.0), 0.75)), 0.25);
    assert_eq!(calc_error_ulp(pred(1.0), near(pred(1.0), 0.75)), 0.75);
    // Infinity is the value after `MAX`
    assert_eq!(calc_error_ulp(f32::MAX, near(f32::MAX, 0.25)), 0.25);
    assert_eq!(calc_error_ulp(f32::INFINITY, near(f32::MAX, 0.25)), 0.75);
    assert_eq!(calc_error_ulp(-f32::INFINITY, near(-f32::MAX, -0.75)), 0.25);
    assert_eq!(calc_error_ulp(f32::INFINITY, f64::INFINITY), 0.0);
    // Subnormals and zero
    let min = f32::from_bits(1);
    assert_eq!(calc_error_ulp(min, 0.0), 1.0);
    assert_eq!(calc_error_ulp(-min, f64::from(min)), 2.0);
    assert_eq!(calc_error_ulp(0.0, f64::from(min) / 4.0), 0.25);
    assert_eq!(
        calc_error_ulp(f32::MIN_POSITIVE, f64::from(pred(f32::MIN_POSITIVE))),
        1.0
    );

    // Exact ties, zeros and NaN
    assert_eq!(calc_error_ulp(1.5, near(1.5, 0.5)), 0.5);
    assert_eq!(calc_error_ulp(succ(1.5), near(1.5, 0.5)), 0.5);
    assert_eq!(calc_error_ulp(0.0, -0.0), 0.0);
    assert_eq!(calc_error_ulp(-0.0, 0.0), 0.0);
    assert_eq!(calc_error_ulp(f32::NAN, f64::NAN), 0.0);
    assert_eq!(calc_error_ulp(1.0, f64::NAN), f32::INFINITY);
    assert_eq!(calc_error_ulp(f32::NAN, 1.0), f32::INFINITY);

    // `expected` rounded to `f32` above the exact result, where the interval
    // below it must be used (the algorithm corrects `lo`): to infinity in
    // `[MAX + ulp/2, 2^128)`, and towards zero to a power of two when
    // negative.
    assert_eq!(calc_error_ulp(f32::MAX, near(f32::MAX, 0.75)), 0.75);
    assert_eq!(calc_error_ulp(f32::INFINITY, near(f32::MAX, 0.75)), 0.25);
    assert_eq!(calc_error_ulp(-1.0, near(-1.0, -0.25)), 0.25);
    assert_eq!(calc_error_ulp(-succ(1.0), near(-1.0, -0.25)), 0.75);
    assert_eq!(calc_error_ulp(-f32::MAX, near(-f32::MAX, -0.25)), 0.25);
    // Tiny negative results, rounded to `-0` (whose ordinal is zero, like
    // `+0`).
    assert_eq!(calc_error_ulp(-0.0, -f64::from(min) / 4.0), 0.25);
    assert_eq!(calc_error_ulp(-min, -f64::from(min) / 4.0), 0.75);
}
