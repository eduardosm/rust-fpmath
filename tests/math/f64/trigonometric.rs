use super::{MIN_MAX_ERROR, RUG_PREC, check_result, values};

#[test]
fn test_sin_cos() {
    test_sin_cos_u(
        Unit::Radians,
        fpmath::sin,
        fpmath::cos,
        fpmath::sin_cos,
        |x| x.sin_cos(rug::Float::new(RUG_PREC)),
    );
}

#[test]
fn test_sind_cosd() {
    test_sin_cos_u(
        Unit::Degrees,
        fpmath::sind,
        fpmath::cosd,
        fpmath::sind_cosd,
        |x| {
            let sin = x.clone().sin_u(360);
            let cos = x.cos_u(360);
            (sin, cos)
        },
    );
}

#[test]
fn test_sinpi_cospi() {
    test_sin_cos_u(
        Unit::HalfRevs,
        fpmath::sinpi,
        fpmath::cospi,
        fpmath::sinpi_cospi,
        |x| {
            let sin = x.clone().sin_pi();
            let cos = x.cos_pi();
            (sin, cos)
        },
    );
}

fn test_sin_cos_u(
    unit: Unit,
    f64_sin: impl Fn(f64) -> f64,
    f64_cos: impl Fn(f64) -> f64,
    f64_sin_cos: impl Fn(f64) -> (f64, f64),
    ref_sin_cos: impl Fn(rug::Float) -> (rug::Float, rug::Float),
) {
    let mut max_sin1_error: f64 = 0.0;
    let mut max_sin2_error: f64 = 0.0;
    let mut max_cos1_error: f64 = 0.0;
    let mut max_cos2_error: f64 = 0.0;
    test_with(unit, |x| {
        let (expected_sin, expected_cos) = ref_sin_cos(rug::Float::with_val(RUG_PREC, x));

        let actual_sin1 = f64_sin(x);
        let actual_cos1 = f64_cos(x);
        let (actual_sin2, actual_cos2) = f64_sin_cos(x);
        assert_result_eq!(f64_sin(-x), -actual_sin1);
        assert_result_eq!(f64_cos(-x), actual_cos1);
        assert_result_eq!(f64_sin_cos(-x), (-actual_sin2, actual_cos2));

        check_result(x, actual_sin1, expected_sin.clone(), &mut max_sin1_error);
        check_result(x, actual_sin2, expected_sin, &mut max_sin2_error);
        check_result(x, actual_cos1, expected_cos.clone(), &mut max_cos1_error);
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
    test_tan_u(Unit::Radians, fpmath::tan, rug::Float::tan);
}

#[test]
fn test_tand() {
    test_tan_u(Unit::Degrees, fpmath::tand, |x| x.tan_u(360));
}

#[test]
fn test_tanpi() {
    test_tan_u(Unit::HalfRevs, fpmath::tanpi, rug::Float::tan_pi);
}

fn test_tan_u(
    unit: Unit,
    f64_tan: impl Fn(f64) -> f64,
    ref_tan: impl Fn(rug::Float) -> rug::Float,
) {
    let mut max_error: f64 = 0.0;
    test_with(unit, |x| {
        let expected = ref_tan(rug::Float::with_val(RUG_PREC, x));
        let actual = f64_tan(x);
        assert_result_eq!(f64_tan(-x), -actual);

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

/// Precision used to calculate the test values.
const EXT_PREC: u32 = 256;

impl Unit {
    /// Returns a right angle (π/2 radians, 90 degrees or 1/2 half revolution)
    /// calculated with precision `prec`.
    fn right_angle(self, prec: u32) -> rug::Float {
        match self {
            Self::Radians => rug::Float::with_val(prec, rug::float::Constant::Pi) / 2,
            Self::Degrees => rug::Float::with_val(prec, 90),
            Self::HalfRevs => rug::Float::with_val(prec, 0.5),
        }
    }

    /// Calculates the sine with high precision.
    fn sin(self, x: &rug::Float) -> rug::Float {
        let x = rug::Float::with_val(EXT_PREC, x);
        match self {
            Self::Radians => x.sin(),
            Self::Degrees => x.sin_u(360),
            Self::HalfRevs => x.sin_pi(),
        }
    }

    /// Calculates the cosine with high precision.
    fn cos(self, x: &rug::Float) -> rug::Float {
        let x = rug::Float::with_val(EXT_PREC, x);
        match self {
            Self::Radians => x.cos(),
            Self::Degrees => x.cos_u(360),
            Self::HalfRevs => x.cos_pi(),
        }
    }

    /// Calculates the arcsine with high precision.
    fn asin(self, y: &rug::Float) -> rug::Float {
        let y = rug::Float::with_val(EXT_PREC, y);
        match self {
            Self::Radians => y.asin(),
            Self::Degrees => y.asin_u(360),
            Self::HalfRevs => y.asin_pi(),
        }
    }

    /// Calculates the arccosine with high precision.
    fn acos(self, y: &rug::Float) -> rug::Float {
        let y = rug::Float::with_val(EXT_PREC, y);
        match self {
            Self::Radians => y.acos(),
            Self::Degrees => y.acos_u(360),
            Self::HalfRevs => y.acos_pi(),
        }
    }
}

fn test_with(unit: Unit, mut f: impl FnMut(f64)) {
    // Negative arguments are tested by the callers.

    // Test special values
    for x in values::specials() {
        f(x);
    }

    // Test a set of subnormal values
    for x in values::subnormals(100_000) {
        f(x);
    }

    // Test across a set of normal values at each binade
    for x in values::binades(-1022..=1023, 1000) {
        f(x);
    }

    // Test values evenly spread over a whole period
    let period = (unit.right_angle(EXT_PREC) * 4u32).to_f64();
    for x in values::interval(0.0, period, 2_000_000) {
        f(x);
    }

    // Test the arguments around the centers (where the reduced argument is
    // zero) and the boundaries (where it is the largest) of the subintervals
    // used by the implementation
    for x in table_points(unit) {
        f(x);
    }

    match unit {
        Unit::Radians => {
            // Test the arguments closest to multiples of a right angle at each
            // binade, which are the hardest to reduce, and whose results are
            // close to zero or infinity
            for x in values::near_multiples(|prec| unit.right_angle(prec), 0..=1023) {
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
    for arg in 1..=200_000 {
        f(arg as f64);
    }

    // Test the arguments around the multiples of small fractions of a right
    // angle (from 1/4 to 1/32), also scaled by 2^20
    for d in [4, 8, 16, 32] {
        let step = unit.right_angle(EXT_PREC) / d;
        for k in 1..=100 {
            for scale in [1, 1 << 20] {
                let x: rug::Float = rug::Float::with_val(EXT_PREC, &step * k) * scale;
                for x in values::around(x.to_f64(), 4) {
                    f(x);
                }
            }
        }
    }
}

/// Returns the arguments around the centers (where the reduced argument is
/// zero) and the boundaries (where it is the largest) of the subintervals
/// used by the implementation, which splits each right angle in 32, in the
/// first two periods and at more binades.
fn table_points(unit: Unit) -> Vec<f64> {
    // Half the width of a subinterval
    let half = unit.right_angle(EXT_PREC) / 64;
    let first = (0..=512).map(|j| rug::Float::with_val(EXT_PREC, j));
    let more =
        values::binades(0..=60, 100).map(|x| rug::Float::with_val(EXT_PREC, x / &half).round());
    first
        .chain(more)
        .flat_map(|j| values::around(rug::Float::with_val(EXT_PREC, j * &half).to_f64(), 2))
        .collect()
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
fn exact_results(unit: Unit) -> Vec<f64> {
    // A step and the number of steps in a period
    let (step, n) = match unit {
        Unit::Radians => return Vec::new(),
        Unit::Degrees => (15.0, 24),
        Unit::HalfRevs => (0.25, 8),
    };
    // The multiples in the first periods, and the multiples in a period
    // scaled by powers of two
    let first = (1..=10_000).map(|j| step * f64::from(j));
    let scaled =
        (1..=n).flat_map(|j| (1..=1023).map(move |k| fpmath::scalbn(step * f64::from(j), k)));
    first.chain(scaled).filter(|x| x.is_finite()).collect()
}

/// Returns arguments whose results are very close to a midpoint between two
/// consecutive values.
///
/// Around their extrema, the sine and the cosine are flat (their results
/// change much less than an ULP between consecutive arguments), so the
/// arguments closest to the inverse of a midpoint have such results (see
/// `values::midpoint_inverse`): small arguments for the cosine and arguments
/// close to a right angle for the sine.
///
/// With small arguments in radians, `sin(x) = x + g(x)` and
/// `tan(x) = x + g(x)` are not flat, but `g` is (see `values::offset_midpoint`).
fn near_midpoints(unit: Unit) -> Vec<f64> {
    let sin = |x: &rug::Float| unit.sin(x);
    let cos = |x: &rug::Float| unit.cos(x);
    let asin = |y: &rug::Float| unit.asin(y);
    let acos = |y: &rug::Float| unit.acos(y);
    let right = unit.right_angle(EXT_PREC);

    // Binades of the distances to the extrema where the functions are flat
    // enough
    let dists = match unit {
        Unit::Radians => -26..=-3,
        Unit::Degrees => -20..=2,
        Unit::HalfRevs => -26..=-5,
    };

    let mut xs = Vec::new();
    for d in values::binades(dists, 500) {
        xs.extend(values::around(values::midpoint_inverse(d, cos, acos), 1));
        let x0 = rug::Float::with_val(EXT_PREC, &right - d).to_f64();
        xs.extend(values::around(values::midpoint_inverse(x0, sin, asin), 1));
    }

    if unit == Unit::Radians {
        let tan = |x: &rug::Float| rug::Float::with_val(EXT_PREC, x).tan();
        for x0 in values::binades(-26..=-1, 500) {
            // g(x) = sin(x) - x, g'(x) = cos(x) - 1
            if let Some(x) = values::offset_midpoint(x0, sin, |x| cos(x) - 1) {
                xs.extend(values::around(x, 1));
            }
            // g(x) = tan(x) - x, g'(x) = tan(x)^2
            if let Some(x) = values::offset_midpoint(x0, tan, |x| tan(x).square()) {
                xs.extend(values::around(x, 1));
            }
        }
    }

    // Close arguments can give the same midpoint
    xs.sort_by(f64::total_cmp);
    xs.dedup();
    xs
}
