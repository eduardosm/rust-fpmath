use super::{MIN_MAX_ERROR, check_result, purify, values};

#[test]
fn test_asin() {
    let mut max_error: f32 = 0.0;
    test_asin_with(Unit::Radians, |x| {
        let expected = fpmath::asin(f64::from(x));
        let actual = fpmath::asin(x);
        assert_result_eq!(fpmath::asin(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_acos() {
    let mut max_error: f32 = 0.0;
    test_acos_with(Unit::Radians, |x| {
        let expected = fpmath::acos(f64::from(x));
        let actual = fpmath::acos(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_asind() {
    let mut max_error: f32 = 0.0;
    test_asin_with(Unit::Degrees, |x| {
        let expected = fpmath::asind(f64::from(x));
        let actual = fpmath::asind(x);
        assert_result_eq!(fpmath::asind(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_acosd() {
    let mut max_error: f32 = 0.0;
    test_acos_with(Unit::Degrees, |x| {
        let expected = fpmath::acosd(f64::from(x));
        let actual = fpmath::acosd(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_asinpi() {
    let mut max_error: f32 = 0.0;
    test_asin_with(Unit::HalfRevs, |x| {
        let expected = fpmath::asinpi(f64::from(x));
        let actual = fpmath::asinpi(x);
        assert_result_eq!(fpmath::asinpi(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_acospi() {
    let mut max_error: f32 = 0.0;
    test_acos_with(Unit::HalfRevs, |x| {
        let expected = fpmath::acospi(f64::from(x));
        let actual = fpmath::acospi(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_atan() {
    let mut max_error: f32 = 0.0;
    test_atan_with(Unit::Radians, |x| {
        let expected = fpmath::atan(f64::from(x));
        let actual = fpmath::atan(x);
        assert_result_eq!(fpmath::atan(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_atand() {
    let mut max_error: f32 = 0.0;
    test_atan_with(Unit::Degrees, |x| {
        let expected = fpmath::atand(f64::from(x));
        let actual = fpmath::atand(x);
        assert_result_eq!(fpmath::atand(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_atanpi() {
    let mut max_error: f32 = 0.0;
    test_atan_with(Unit::HalfRevs, |x| {
        let expected = fpmath::atanpi(f64::from(x));
        let actual = fpmath::atanpi(x);
        assert_result_eq!(fpmath::atanpi(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_atan2() {
    let mut max_error: f32 = 0.0;
    test_atan2_with(|y, x| {
        let expected = fpmath::atan2(f64::from(y), f64::from(x));
        let actual = fpmath::atan2(y, x);
        assert_result_eq!(fpmath::atan2(-y, x), -actual);

        check_result((y, x), actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_atan2d() {
    let mut max_error: f32 = 0.0;
    test_atan2_with(|y, x| {
        let expected = fpmath::atan2d(f64::from(y), f64::from(x));
        let actual = fpmath::atan2d(y, x);
        assert_result_eq!(fpmath::atan2d(-y, x), -actual);

        check_result((y, x), actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_atan2pi() {
    let mut max_error: f32 = 0.0;
    test_atan2_with(|y, x| {
        let expected = fpmath::atan2pi(f64::from(y), f64::from(x));
        let actual = fpmath::atan2pi(y, x);
        assert_result_eq!(fpmath::atan2pi(-y, x), -actual);

        check_result((y, x), actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

/// Angle unit of the results of the inverse trigonometric functions.
#[derive(Copy, Clone, PartialEq, Eq)]
enum Unit {
    Radians,
    Degrees,
    HalfRevs,
}

impl Unit {
    /// Calculates the cosine in `f64`.
    fn cos(self, x: f64) -> f64 {
        match self {
            Self::Radians => fpmath::cos(x),
            Self::Degrees => fpmath::cosd(x),
            Self::HalfRevs => fpmath::cospi(x),
        }
    }

    /// Calculates the arccosine in `f64`.
    fn acos(self, x: f64) -> f64 {
        match self {
            Self::Radians => fpmath::acos(x),
            Self::Degrees => fpmath::acosd(x),
            Self::HalfRevs => fpmath::acospi(x),
        }
    }

    /// Calculates the tangent in `f64`.
    fn tan(self, x: f64) -> f64 {
        match self {
            Self::Radians => fpmath::tan(x),
            Self::Degrees => fpmath::tand(x),
            Self::HalfRevs => fpmath::tanpi(x),
        }
    }

    /// Calculates the arctangent in `f64`.
    fn atan(self, x: f64) -> f64 {
        match self {
            Self::Radians => fpmath::atan(x),
            Self::Degrees => fpmath::atand(x),
            Self::HalfRevs => fpmath::atanpi(x),
        }
    }
}

/// Generates the arguments for `asin`, and the non-negative arguments for
/// `acos`.
fn asin_args(mut f: impl FnMut(f32)) {
    // Test special values
    for x in values::specials() {
        f(x);
    }

    // Test a set of subnormal numbers (including all the bit lengths), whose
    // results are trivial (`x`, a right angle, or proportional to `x`)
    for x in values::subnormals(10_000) {
        f(x);
    }

    // Test across a wide range of normal numbers below one, and the mantissa
    // patterns at the binades above one (whose results are NaN)
    for x in values::binades(-126..=-1, 10_000).chain(values::binades(0..=127, 0)) {
        f(x);
    }

    // Test more values at the binades where the results are not just
    // proportional to `x`
    for x in values::binades(-14..=-2, 200_000) {
        f(x);
    }

    // Exhaustive test of all the values in [1/2, 1], where the implementation
    // uses a different evaluation method
    for x in values::interval_full(0.5, 1.0) {
        f(x);
    }

    for i in 1..=1000 {
        let x = purify((i as f32) / 1000.0);
        f(x);
    }
}

fn test_asin_with(unit: Unit, mut f: impl FnMut(f32)) {
    // Negative arguments are tested by the callers.

    asin_args(&mut f);

    if unit == Unit::Radians {
        // Test results that are very close to a midpoint between two
        // consecutive values, which are the hardest to round. With small
        // arguments, `asin(x) = x + g(x)` is not flat, but `g` is (see
        // `values::offset_midpoint`).
        // g(x) = asin(x) - x, g'(x) = 1 / sqrt(1 - x^2) - 1
        let dg = |x: f64| 1.0 / (1.0 - x * x).sqrt() - 1.0;
        for x0 in values::binades(-12..=-2, 20_000) {
            if let Some(x) = values::offset_midpoint(x0, fpmath::asin, dg) {
                for x in values::around(x, 1) {
                    f(x);
                }
            }
        }
    }
}

fn test_acos_with(unit: Unit, mut f: impl FnMut(f32)) {
    asin_args(|x| {
        f(x);
        f(-x);
    });

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round. With small arguments, the
    // results are close to a right angle and the function is flat, so the
    // arguments closest to the inverse of a midpoint have such results (see
    // `values::midpoint_inverse`).
    let acos = |x: f64| unit.acos(x);
    let cos = |y: f64| unit.cos(y);
    for x0 in values::binades(-20..=-3, 20_000) {
        for x0 in [x0, -x0] {
            for x in values::around(values::midpoint_inverse(x0, acos, cos), 1) {
                f(x);
            }
        }
    }
}

fn test_atan_with(unit: Unit, mut f: impl FnMut(f32)) {
    // Negative arguments are tested by the callers.

    // Test special values
    for x in values::specials() {
        f(x);
    }

    // Test a set of subnormal numbers (including all the bit lengths), whose
    // results are trivial (`x`, a right angle, or proportional to `x`)
    for x in values::subnormals(10_000) {
        f(x);
    }

    // Test across a wide range of normal numbers
    for x in values::binades(-126..=127, 10_000) {
        f(x);
    }

    // Test more values at the binades where the results are not just
    // proportional to `x` or a right angle
    for x in values::binades(-14..=26, 100_000) {
        f(x);
    }

    // Exhaustive test of all the values in [1/2, 2), around one, where the
    // implementation changes the evaluation method
    for x in values::binades_full(-1..=0) {
        f(x);
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round. With large arguments, the
    // results are close to a right angle and the function is flat, so the
    // arguments closest to the inverse of a midpoint have such results (see
    // `values::midpoint_inverse`).
    let atan = |x: f64| unit.atan(x);
    let tan = |y: f64| unit.tan(y);
    for x0 in values::binades(3..=26, 20_000) {
        for x in values::around(values::midpoint_inverse(x0, atan, tan), 1) {
            f(x);
        }
    }

    if unit == Unit::Radians {
        // With small arguments, `atan(x) = x + g(x)` is not flat, but `g` is
        // (see `values::offset_midpoint`).
        // g(x) = atan(x) - x, g'(x) = 1 / (1 + x^2) - 1 = -x^2 / (1 + x^2)
        let dg = |x: f64| -(x * x) / (1.0 + x * x);
        for x0 in values::binades(-12..=-1, 20_000) {
            if let Some(x) = values::offset_midpoint(x0, fpmath::atan, dg) {
                for x in values::around(x, 1) {
                    f(x);
                }
            }
        }
    }
}

fn test_atan2_with(mut f: impl FnMut(f32, f32)) {
    // Toggling the sign of `y` is done by the callers.

    // Test special values
    for y in values::specials() {
        for x in values::specials() {
            f(y, x);
            f(y, -x);
        }
    }

    // Test the ratios `|y / x|` between 2^-64 and 2^64, which include all
    // the ones where `atan2(y, x)` is not just `y / x` (below 2^-32) or `±π/2`
    // (above 2^24), at every scale, including subnormals
    for (y, x) in values::binade_pairs(values::exp_pairs(-64..=64, 1), 20) {
        f(y, x);
        f(y, -x);
    }

    // Test more values with ratios between 2^-8 and 2^8, where the
    // approximation of `atan` is less accurate, at a few scales (since only the
    // ratio matters away from the limits of the range)
    for (y, x) in values::binade_pairs(values::exp_pairs(-8..=8, 16), 20_000) {
        f(y, x);
        f(y, -x);
    }

    // Test larger ratios at many scales
    let large_diffs = (65..=276).flat_map(|d| [d, -d]);
    for (y, x) in values::binade_pairs(values::exp_pairs(large_diffs, 4), 1) {
        f(y, x);
        f(y, -x);
    }

    // Test a set of subnormal numbers (including all the bit lengths) with one
    // as the other argument, whose results are trivial
    for v in values::subnormals(10_000) {
        f(v, 1.0);
        f(v, -1.0);
        f(1.0, v);
        f(1.0, -v);
    }

    // Test `|y|` close to `|x|`, where the arguments are swapped
    for x in values::binades(-126..=127, 1000).chain(values::subnormals(10_000)) {
        for y in values::around(x, 1) {
            f(y, x);
            f(y, -x);
        }
    }
}
