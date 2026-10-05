use super::{MIN_MAX_ERROR, check_result, values};

#[test]
fn test_sin_cos() {
    test_sin_cos_u(
        Unit::Radians,
        fpmath::sin,
        fpmath::cos,
        fpmath::sin_cos,
        fpmath::sin_cos,
    );
}

#[test]
fn test_sind_cosd() {
    test_sin_cos_u(
        Unit::Degrees,
        fpmath::sind,
        fpmath::cosd,
        fpmath::sind_cosd,
        fpmath::sind_cosd,
    );
}

#[test]
fn test_sinpi_cospi() {
    test_sin_cos_u(
        Unit::HalfRevs,
        fpmath::sinpi,
        fpmath::cospi,
        fpmath::sinpi_cospi,
        fpmath::sinpi_cospi,
    );
}

fn test_sin_cos_u(
    unit: Unit,
    f32_sin: impl Fn(f32) -> f32,
    f32_cos: impl Fn(f32) -> f32,
    f32_sin_cos: impl Fn(f32) -> (f32, f32),
    ref_sin_cos: impl Fn(f64) -> (f64, f64),
) {
    let mut max_sin1_error: f32 = 0.0;
    let mut max_sin2_error: f32 = 0.0;
    let mut max_cos1_error: f32 = 0.0;
    let mut max_cos2_error: f32 = 0.0;
    test_with(unit, |x| {
        let (expected_sin, expected_cos) = ref_sin_cos(f64::from(x));

        let actual_sin1 = f32_sin(x);
        let actual_cos1 = f32_cos(x);
        let (actual_sin2, actual_cos2) = f32_sin_cos(x);
        assert_result_eq!(f32_sin(-x), -actual_sin1);
        assert_result_eq!(f32_cos(-x), actual_cos1);
        assert_result_eq!(f32_sin_cos(-x), (-actual_sin2, actual_cos2));

        check_result(x, actual_sin1, expected_sin, &mut max_sin1_error);
        check_result(x, actual_sin2, expected_sin, &mut max_sin2_error);
        check_result(x, actual_cos1, expected_cos, &mut max_cos1_error);
        check_result(x, actual_cos2, expected_cos, &mut max_cos2_error);
    });
    eprintln!("max sin1 error = {max_sin1_error}");
    eprintln!("max sin2 error = {max_sin2_error}");
    eprintln!("max cos1 error = {max_cos1_error}");
    eprintln!("max cos2 error = {max_cos2_error}");
    assert!(max_sin1_error > MIN_MAX_ERROR);
    assert!(max_sin2_error > MIN_MAX_ERROR);
    assert!(max_cos1_error > MIN_MAX_ERROR);
    assert!(max_cos2_error > MIN_MAX_ERROR);
}

#[test]
fn test_tan() {
    test_tan_u(Unit::Radians, fpmath::tan, fpmath::tan);
}

#[test]
fn test_tand() {
    test_tan_u(Unit::Degrees, fpmath::tand, fpmath::tand);
}

#[test]
fn test_tanpi() {
    test_tan_u(Unit::HalfRevs, fpmath::tanpi, fpmath::tanpi);
}

fn test_tan_u(unit: Unit, f32_tan: impl Fn(f32) -> f32, ref_tan: impl Fn(f64) -> f64) {
    let mut max_error: f32 = 0.0;
    test_with(unit, |x| {
        let expected = ref_tan(f64::from(x));
        let actual = f32_tan(x);
        assert_result_eq!(f32_tan(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

/// Angle unit of the trigonometric functions.
#[derive(Copy, Clone, PartialEq, Eq)]
enum Unit {
    Radians,
    Degrees,
    HalfRevs,
}

impl Unit {
    /// Returns a right angle (π/2 radians, 90 degrees or 1/2 half revolution).
    fn right_angle(self) -> f64 {
        match self {
            Self::Radians => std::f64::consts::FRAC_PI_2,
            Self::Degrees => 90.0,
            Self::HalfRevs => 0.5,
        }
    }

    /// Calculates the cosine in `f64`.
    fn cos(self, x: f64) -> f64 {
        match self {
            Self::Radians => fpmath::cos(x),
            Self::Degrees => fpmath::cosd(x),
            Self::HalfRevs => fpmath::cospi(x),
        }
    }

    /// Calculates the arccosine in `f64`.
    fn acos(self, y: f64) -> f64 {
        match self {
            Self::Radians => fpmath::acos(y),
            Self::Degrees => fpmath::acosd(y),
            Self::HalfRevs => fpmath::acospi(y),
        }
    }
}

fn test_with(unit: Unit, mut f: impl FnMut(f32)) {
    // Negative arguments are tested by the callers.

    // Test special values
    for x in values::specials() {
        f(x);
    }

    // Exhaustive test of all subnormal numbers
    for x in values::subnormals_full() {
        f(x);
    }

    // Test across a wide range of normal numbers
    for x in values::binades(-126..=127, 20_000) {
        f(x);
    }

    // Exhaustive test of all the values in a whole period (except the small
    // ones)
    let period = match unit {
        // [1/2, 8)
        Unit::Radians => -1..=2,
        // [32, 512)
        Unit::Degrees => 5..=8,
        // [1/4, 2)
        Unit::HalfRevs => -2..=0,
    };
    for x in values::binades_full(period) {
        f(x);
    }

    // Test the arguments around the centers (where the reduced argument is
    // zero) and the boundaries (where it is the largest) of the subintervals
    // used by the implementation, beyond the first period
    for x in table_points(unit) {
        f(x);
    }

    match unit {
        Unit::Radians => {
            // Test the arguments closest to multiples of a right angle at each
            // binade, which are the hardest to reduce, and whose results are
            // close to zero or infinity
            let right_angle = |prec| rug::Float::with_val(prec, rug::float::Constant::Pi) / 2;
            for x in values::near_multiples(right_angle, 0..=127) {
                f(x);
            }
        }
        Unit::Degrees | Unit::HalfRevs => {
            // Test the arguments with exact results, at all magnitudes
            for x in exact_results(unit) {
                f(x);
            }
        }
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round
    for x in near_midpoints(unit) {
        f(x);
    }

    // Test integers
    for arg in 1..=20_000 {
        f(arg as f32);
    }

    // Test the arguments around the multiples of small fractions of a right
    // angle (from 1/4 to 1/32), also scaled by 2^9
    for d in [4.0, 8.0, 16.0, 32.0] {
        let step = unit.right_angle() / d;
        for k in 1..=100 {
            for scale in [1.0, f64::from(1 << 9)] {
                let x = (step * f64::from(k) * scale) as f32;
                for x in values::around(x, 4) {
                    f(x);
                }
            }
        }
    }
}

/// Returns the arguments around the centers (where the reduced argument is
/// zero) and the boundaries (where it is the largest) of the subintervals
/// used by the implementation, which splits each right angle in 32, at many
/// binades.
fn table_points(unit: Unit) -> impl Iterator<Item = f32> {
    // Half the width of a subinterval
    let half = unit.right_angle() / 64.0;
    values::binades(0..=40, 1000).flat_map(move |x| {
        let j = (f64::from(x) / half).round();
        values::around((j * half) as f32, 2)
    })
}

/// Returns arguments that include the ones with exact results, at all
/// magnitudes: multiples of 15 degrees and of a quarter of half revolution.
///
/// In degrees, the sine is exact (`0`, `±1/2` or `±1`) at the multiples of 30
/// degrees other than `±60` and `±120` modulo 360, the cosine at the multiples
/// of 30 degrees other than `±30` and `±150` modulo 360, and the tangent (`0`,
/// `±1` or infinite) at the multiples of 45 degrees. In half revolutions, the
/// sine and the cosine are exact (`0` or `±1`) at the multiples of a half, and
/// the tangent at the multiples of a quarter.
fn exact_results(unit: Unit) -> Vec<f32> {
    // A step and the number of steps in a period
    let (step, n) = match unit {
        Unit::Radians => return Vec::new(),
        Unit::Degrees => (15.0, 24),
        Unit::HalfRevs => (0.25, 8),
    };
    // The multiples in the first periods, and the multiples in a period
    // scaled by powers of two
    let first = (1..=10_000).map(|j| step * j as f32);
    let scaled = (1..=n).flat_map(|j| (1..=127).map(move |k| fpmath::scalbn(step * j as f32, k)));
    first.chain(scaled).filter(|x| x.is_finite()).collect()
}

/// Returns arguments whose results are very close to a midpoint between two
/// consecutive values.
///
/// Around zero, the cosine is flat (its result changes much less than an ULP
/// between consecutive arguments), so the arguments closest to the inverse of
/// a midpoint have such results (see `values::midpoint_inverse`). Other
/// extrema of the sine and the cosine are in the period that is tested
/// exhaustively.
///
/// With small arguments in radians, `sin(x) = x + g(x)` and
/// `tan(x) = x + g(x)` are not flat, but `g` is (see `values::offset_midpoint`).
fn near_midpoints(unit: Unit) -> Vec<f32> {
    // Binades of the arguments where the cosine is flat enough
    let small = match unit {
        Unit::Radians => -13..=-3,
        Unit::Degrees => -7..=2,
        Unit::HalfRevs => -13..=-5,
    };

    let mut xs = Vec::new();
    for x0 in values::binades(small, 20_000) {
        let x = values::midpoint_inverse(x0, |x| unit.cos(x), |y| unit.acos(y));
        xs.extend(values::around(x, 1));
    }

    if unit == Unit::Radians {
        for x0 in values::binades(-11..=-1, 20_000) {
            // g(x) = sin(x) - x, g'(x) = cos(x) - 1
            if let Some(x) = values::offset_midpoint(x0, fpmath::sin, |x| fpmath::cos(x) - 1.0) {
                xs.extend(values::around(x, 1));
            }
            // g(x) = tan(x) - x, g'(x) = tan(x)^2
            let dg = |x: f64| fpmath::tan(x) * fpmath::tan(x);
            if let Some(x) = values::offset_midpoint(x0, fpmath::tan, dg) {
                xs.extend(values::around(x, 1));
            }
        }
    }

    // Close arguments can give the same midpoint
    xs.sort_by(f32::total_cmp);
    xs.dedup();
    xs
}
