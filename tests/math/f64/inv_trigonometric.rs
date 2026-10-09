use super::{MIN_MAX_ERROR, RUG_PREC, check_result, purify, values};

#[test]
fn test_asin() {
    let mut max_error: f64 = 0.0;
    test_asin_with(Unit::Radians, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).asin();
        let actual = fpmath::asin(x);
        assert_result_eq!(fpmath::asin(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_acos() {
    let mut max_error: f64 = 0.0;
    test_acos_with(Unit::Radians, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).acos();
        let actual = fpmath::acos(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_asind() {
    let mut max_error: f64 = 0.0;
    test_asin_with(Unit::Degrees, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).asin_u(360);
        let actual = fpmath::asind(x);
        assert_result_eq!(fpmath::asind(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_acosd() {
    let mut max_error: f64 = 0.0;
    test_acos_with(Unit::Degrees, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).acos_u(360);
        let actual = fpmath::acosd(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_asinpi() {
    let mut max_error: f64 = 0.0;
    test_asin_with(Unit::HalfRevs, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).asin_pi();
        let actual = fpmath::asinpi(x);
        assert_result_eq!(fpmath::asinpi(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_acospi() {
    let mut max_error: f64 = 0.0;
    test_acos_with(Unit::HalfRevs, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).acos_pi();
        let actual = fpmath::acospi(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_atan() {
    let mut max_error: f64 = 0.0;
    test_atan_with(Unit::Radians, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).atan();
        let actual = fpmath::atan(x);
        assert_result_eq!(fpmath::atan(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_atand() {
    let mut max_error: f64 = 0.0;
    test_atan_with(Unit::Degrees, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).atan_u(360);
        let actual = fpmath::atand(x);
        assert_result_eq!(fpmath::atand(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_atanpi() {
    let mut max_error: f64 = 0.0;
    test_atan_with(Unit::HalfRevs, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).atan_pi();
        let actual = fpmath::atanpi(x);
        assert_result_eq!(fpmath::atanpi(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_atan2() {
    let mut max_error: f64 = 0.0;
    test_atan2_with(|y, x| {
        let expected = rug::Float::with_val(RUG_PREC, y).atan2(&rug::Float::with_val(RUG_PREC, x));
        let actual = fpmath::atan2(y, x);
        assert_result_eq!(fpmath::atan2(-y, x), -actual);

        check_result((y, x), actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_atan2d() {
    let mut max_error: f64 = 0.0;
    test_atan2_with(|y, x| {
        let expected =
            rug::Float::with_val(RUG_PREC, y).atan2_u(&rug::Float::with_val(RUG_PREC, x), 360);
        let actual = fpmath::atan2d(y, x);
        assert_result_eq!(fpmath::atan2d(-y, x), -actual);

        check_result((y, x), actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_atan2pi() {
    let mut max_error: f64 = 0.0;
    test_atan2_with(|y, x| {
        let expected =
            rug::Float::with_val(RUG_PREC, y).atan2_pi(&rug::Float::with_val(RUG_PREC, x));
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

/// Precision used to calculate the test values.
const EXT_PREC: u32 = 256;

impl Unit {
    /// Calculates the cosine with high precision.
    fn cos(self, x: &rug::Float) -> rug::Float {
        let x = rug::Float::with_val(EXT_PREC, x);
        match self {
            Self::Radians => x.cos(),
            Self::Degrees => x.cos_u(360),
            Self::HalfRevs => x.cos_pi(),
        }
    }

    /// Calculates the arccosine with high precision.
    fn acos(self, x: &rug::Float) -> rug::Float {
        let x = rug::Float::with_val(EXT_PREC, x);
        match self {
            Self::Radians => x.acos(),
            Self::Degrees => x.acos_u(360),
            Self::HalfRevs => x.acos_pi(),
        }
    }

    /// Calculates the tangent with high precision.
    fn tan(self, x: &rug::Float) -> rug::Float {
        let x = rug::Float::with_val(EXT_PREC, x);
        match self {
            Self::Radians => x.tan(),
            Self::Degrees => x.tan_u(360),
            Self::HalfRevs => x.tan_pi(),
        }
    }

    /// Calculates the arctangent with high precision.
    fn atan(self, x: &rug::Float) -> rug::Float {
        let x = rug::Float::with_val(EXT_PREC, x);
        match self {
            Self::Radians => x.atan(),
            Self::Degrees => x.atan_u(360),
            Self::HalfRevs => x.atan_pi(),
        }
    }
}

/// Generates the arguments for `asin`, and the non-negative arguments for
/// `acos`.
fn asin_args(mut f: impl FnMut(f64)) {
    // Test special values
    for x in values::specials() {
        f(x);
    }

    // Test a set of subnormal values
    for x in values::subnormals(100_000) {
        f(x);
    }

    // Test a set of normal values at each binade below one, and the mantissa
    // patterns at the binades above one (whose results are NaN)
    for x in values::binades(-1022..=-1, 1000).chain(values::binades(0..=1023, 0)) {
        f(x);
    }

    // Test more values at the binades where the results are not just
    // proportional to `x`
    for x in values::binades(-30..=-1, 50_000) {
        f(x);
    }

    // Test all the values within 2^16 ULPs of one and one half, where the
    // implementation changes the evaluation method
    for x in values::around(1.0, 1 << 16).chain(values::around(0.5, 1 << 16)) {
        f(x);
    }

    // Test the values around 2^-32, where the implementation changes the
    // evaluation method (for tiny arguments)
    for x in values::around(fpmath::scalbn(1.0, -32), 10_000) {
        f(x);
    }

    // Test values close to one at larger distances (`1 - v` with `v` at each
    // binade below one half)
    for v in values::binades(-37..=-2, 10_000) {
        f(purify(1.0 - v));
    }

    for i in 1..=1000 {
        let x = purify((i as f64) / 1000.0);
        f(x);
    }
}

fn test_asin_with(unit: Unit, mut f: impl FnMut(f64)) {
    // Negative arguments are tested by the callers.

    asin_args(&mut f);

    if unit == Unit::Radians {
        // Test results that are very close to a midpoint between two
        // consecutive values, which are the hardest to round. With small
        // arguments, `asin(x) = x + g(x)` is not flat, but `g` is (see
        // `values::offset_midpoint`).
        let asin = |x: &rug::Float| rug::Float::with_val(EXT_PREC, x).asin();
        // g(x) = asin(x) - x, g'(x) = 1 / sqrt(1 - x^2) - 1
        let dg = |x: &rug::Float| -> rug::Float {
            let r = rug::Float::with_val(EXT_PREC, 1 - x.clone().square());
            r.sqrt().recip() - 1
        };
        for x0 in values::binades(-26..=-2, 1000) {
            if let Some(x) = values::offset_midpoint(x0, asin, dg) {
                for x in values::around(x, 1) {
                    f(x);
                }
            }
        }
    }
}

fn test_acos_with(unit: Unit, mut f: impl FnMut(f64)) {
    asin_args(|x| {
        f(x);
        f(-x);
    });

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round. With small arguments, the
    // results are close to a right angle and the function is flat, so the
    // arguments closest to the inverse of a midpoint have such results (see
    // `values::midpoint_inverse`).
    let acos = |x: &rug::Float| unit.acos(x);
    let cos = |y: &rug::Float| unit.cos(y);
    for x0 in values::binades(-40..=-3, 1000) {
        for x0 in [x0, -x0] {
            for x in values::around(values::midpoint_inverse(x0, acos, cos), 1) {
                f(x);
            }
        }
    }
}

fn test_atan_with(unit: Unit, mut f: impl FnMut(f64)) {
    // Negative arguments are tested by the callers.

    // Test special values
    for x in values::specials() {
        f(x);
    }

    // Test a set of subnormal values
    for x in values::subnormals(100_000) {
        f(x);
    }

    // Test a set of normal values at each binade
    for x in values::binades(-1022..=1023, 1000) {
        f(x);
    }

    // Test more values at the binades where the results are not just
    // proportional to `x` or a right angle
    for x in values::binades(-30..=55, 20_000) {
        f(x);
    }

    // Test all the values within 2^16 ULPs of one, where the implementation
    // changes the evaluation method
    for x in values::around(1.0, 1 << 16) {
        f(x);
    }

    // Test the values around 2^-32 and 2^60, where the implementation also
    // changes the evaluation method (for tiny and huge arguments)
    for x in values::around(fpmath::scalbn(1.0, -32), 10_000)
        .chain(values::around(fpmath::scalbn(1.0, 60), 10_000))
    {
        f(x);
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round. With large arguments, the
    // results are close to a right angle and the function is flat, so the
    // arguments closest to the inverse of a midpoint have such results (see
    // `values::midpoint_inverse`).
    let atan = |x: &rug::Float| unit.atan(x);
    let tan = |y: &rug::Float| unit.tan(y);
    for x0 in values::binades(3..=55, 1000) {
        for x in values::around(values::midpoint_inverse(x0, atan, tan), 1) {
            f(x);
        }
    }

    if unit == Unit::Radians {
        // With small arguments, `atan(x) = x + g(x)` is not flat, but `g` is
        // (see `values::offset_midpoint`).
        // g(x) = atan(x) - x, g'(x) = 1 / (1 + x^2) - 1 = -x^2 / (1 + x^2)
        let dg = |x: &rug::Float| -> rug::Float {
            let x2 = rug::Float::with_val(EXT_PREC, x.square_ref());
            -(rug::Float::with_val(EXT_PREC, &x2 / (x2.clone() + 1)))
        };
        for x0 in values::binades(-26..=-1, 1000) {
            if let Some(x) = values::offset_midpoint(x0, atan, dg) {
                for x in values::around(x, 1) {
                    f(x);
                }
            }
        }
    }
}

fn test_atan2_with(mut f: impl FnMut(f64, f64)) {
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
    // (above 2^53), at many scales, including subnormals
    for (y, x) in values::binade_pairs(values::exp_pairs(-64..=64, 32), 4) {
        f(y, x);
        f(y, -x);
    }

    // Test more values with ratios between 2^-8 and 2^8, where the
    // approximation of `atan` is less accurate, at a few scales (since only the
    // ratio matters away from the limits of the range)
    for (y, x) in values::binade_pairs(values::exp_pairs(-8..=8, 256), 20_000) {
        f(y, x);
        f(y, -x);
    }

    // Test larger ratios at fewer scales, including the ones where the
    // quotient is subnormal and where scaling a subnormal argument makes the
    // other one overflow
    let large_diffs = (65..=2097).step_by(18).flat_map(|d| [d, -d]);
    for (y, x) in values::binade_pairs(values::exp_pairs(large_diffs, 61), 1) {
        f(y, x);
        f(y, -x);
    }

    // Test the larger argument in the binades around 2^1021, where the
    // implementation starts scaling the arguments (to avoid an overflow and a
    // subnormal reciprocal), with ratios `|y / x|` from 2^-33 to 2^33 (which
    // include the limits of the table-based reduction at 2^-32 and 2^32)
    let near_scaling: Vec<(i16, i16)> = (1019..=1023)
        .flat_map(|e| {
            [0, 1, 2, 6, 7, 31, 32]
                .into_iter()
                .flat_map(move |d| [(e, e - d), (e - d, e)])
        })
        .collect();
    for (y, x) in values::binade_pairs(near_scaling, 100) {
        f(y, x);
        f(y, -x);
    }

    // Test subnormal numbers with one as the other argument
    for v in values::subnormals(10_000) {
        f(v, 1.0);
        f(v, -1.0);
        f(1.0, v);
        f(1.0, -v);
    }

    // Test `|y|` close to `|x|`, where the arguments are swapped
    for x in values::binades(-1022..=1023, 8).chain(values::subnormals(1000)) {
        for y in values::around(x, 1) {
            f(y, x);
            f(y, -x);
        }
    }
}
