use super::{MIN_MAX_ERROR, check_result, values};

#[test]
fn test_erf() {
    let mut max_error: f32 = 0.0;
    test_with(Func::Erf, |x| {
        let expected = fpmath::erf(f64::from(x));
        let actual = fpmath::erf(x);
        assert_result_eq!(fpmath::erf(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_erfc() {
    let mut max_error: f32 = 0.0;
    test_with(Func::Erfc, |x| {
        let expected = fpmath::erfc(f64::from(x));
        let actual = fpmath::erfc(x);
        check_result(x, actual, expected, &mut max_error);

        let expected = fpmath::erfc(f64::from(-x));
        let actual = fpmath::erfc(-x);
        check_result(-x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_erfcx() {
    let mut max_error: f32 = 0.0;
    test_erfcx_with(|x| {
        let expected = fpmath::erfcx(f64::from(x));
        let actual = fpmath::erfcx(x);
        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

/// An error function.
#[derive(Copy, Clone, PartialEq, Eq)]
enum Func {
    Erf,
    Erfc,
}

impl Func {
    /// Calculates the function in `f64`.
    fn eval(self, x: f64) -> f64 {
        match self {
            Self::Erf => fpmath::erf(x),
            Self::Erfc => fpmath::erfc(x),
        }
    }

    /// Solves `func(x) = y` in `f64` with Newton's method, starting from
    /// `x0`.
    fn solve(self, y: f64, x0: f64) -> f64 {
        // erf'(x) = 2 / sqrt(π) * exp(-x^2) and erfc'(x) = -erf'(x)
        let sign = match self {
            Self::Erf => 1.0,
            Self::Erfc => -1.0,
        };
        let mut x = x0;
        for _ in 0..100 {
            let d = sign * std::f64::consts::FRAC_2_SQRT_PI * fpmath::exp(-x * x);
            let dx = (self.eval(x) - y) / d;
            x -= dx;
            if dx.abs() <= x.abs() * f64::EPSILON {
                break;
            }
        }
        x
    }
}

/// Arguments from which `erf` rounds to one and `erfc` (of the negated
/// argument) rounds to two.
const SATURATION_MAX: f32 = 4.0;

/// Arguments from which `erfc` rounds to zero.
const ERFC_ZERO_MAX: f32 = 10.25;

fn test_with(func: Func, mut f: impl FnMut(f32)) {
    // Negative arguments are tested by the callers.

    // The largest argument with a non-trivial result
    let x_max = match func {
        Func::Erf => SATURATION_MAX,
        Func::Erfc => ERFC_ZERO_MAX,
    };

    // Test special values
    for x in values::specials() {
        f(x);
    }

    // Test a set of subnormal values, whose results are subnormal with `erf`
    for x in values::subnormals(100_000) {
        f(x);
    }

    // Test a set of normal values at each binade. The results of the large
    // ones are trivial (zero, one or two), and so are the results of `erfc`
    // with the tiny ones (one).
    for x in values::binades(-126..=127, 1000) {
        f(x);
    }

    // Test more values with small magnitudes, whose results are calculated
    // with the same polynomial (only with its first term below 2^-27, the
    // values closest to each power of two are included)
    for x in values::binades(-30..=-2, 200_000) {
        f(x);
    }

    // Test values evenly spread over the whole range of non-trivial results
    for x in values::interval(0.0, x_max + 0.25, 20_000_000) {
        f(x);
    }

    // Exhaustive test of all the arguments from 0.5 up to the saturation of
    // `erf` and of `erfc` with negative arguments
    for x in values::binades_full(-1..=1) {
        f(x);
    }

    // Exhaustive test of all the arguments with subnormal results
    if func == Func::Erfc {
        let (lo, hi) = erfc_subnormal_range();
        for x in values::interval_full(lo, hi) {
            f(x);
        }
    }

    // Test the arguments around the limits of the range
    for x in limits() {
        for x in values::around(x, 100_000) {
            f(x);
        }
    }

    // Test the arguments around the limits of the intervals of the
    // approximations of `erfcx`
    for x in interval_limits() {
        for x in values::around(x, 100) {
            f(x);
        }
    }

    // Test the arguments around the boundaries of the argument reduction of
    // the exponential
    for x in reduction_boundaries(x_max) {
        for x in values::around(x, 2) {
            f(x);
        }
    }

    // Test results of `erfc` that are very close to a midpoint between two
    // consecutive values, which are the hardest to round
    if func == Func::Erfc {
        for x in erfc_near_midpoints() {
            f(x);
        }
    }
}

/// Returns the range of arguments where the results of `erfc` are subnormal
/// and not zero.
fn erfc_subnormal_range() -> (f32, f32) {
    (
        Func::Erfc.solve(fpmath::scalbn(1.0, -126), 9.1) as f32,
        Func::Erfc.solve(fpmath::scalbn(1.0, -150), 10.0) as f32,
    )
}

/// Returns the arguments at the limits of the range: where the functions
/// change the evaluation method (0.5 and the saturation limits), where `erfc`
/// rounds to one with tiny arguments (positive and negative), where `erf`
/// rounds to one and `erfc` rounds to two (with the negated argument), and
/// where `erfc` becomes subnormal and rounds to zero.
fn limits() -> [f32; 9] {
    let exp2i = |e: i32| fpmath::scalbn(1.0f64, e);
    let erfc_inv = |e: i32, x0: f64| Func::Erfc.solve(exp2i(e), x0) as f32;
    // erf(x) ~= x * 2 / sqrt(π) with tiny x
    let tiny = |e: i32| (std::f64::consts::PI.sqrt() * exp2i(e - 1)) as f32;
    let (subnormal, zero) = erfc_subnormal_range();
    [
        0.5,
        SATURATION_MAX,
        ERFC_ZERO_MAX,
        tiny(-25),
        tiny(-24),
        erfc_inv(-25, 3.9),
        erfc_inv(-24, 3.8),
        subnormal,
        zero,
    ]
}

/// Returns the limits of the intervals of the approximations of `erfcx`:
/// each binade from 0.5 to 16 is split in four intervals of the same length.
fn interval_limits() -> Vec<f32> {
    (-1..=3)
        .flat_map(|e| (0..4).map(move |j| fpmath::scalbn(1.0 + 0.25 * j as f32, e)))
        .collect()
}

/// Number of parts in which the argument reduction of `exp` splits each power
/// of two: the arguments are reduced to `m * ln(2) / N + r`, with an integer
/// `m` and `|r| <= ln(2) / (2 * N)`.
const REDUCTION_N: u32 = 128;

/// Returns the arguments up to `x_max` where `m` changes in the argument
/// reduction of `exp(-x^2)`, `sqrt((m + 1/2) * ln(2) / N)`.
fn reduction_boundaries(x_max: f32) -> impl Iterator<Item = f32> {
    let x_max = f64::from(x_max);
    let m_max = (x_max * x_max / std::f64::consts::LN_2 * f64::from(REDUCTION_N)) as u32;
    (0..=m_max).map(|m| {
        let x2 = f64::from(2 * m + 1) * std::f64::consts::LN_2 / f64::from(2 * REDUCTION_N);
        x2.sqrt() as f32
    })
}

/// Returns arguments whose results with `erfc` (and with `erfc` of the negated
/// argument) are very close to a midpoint between two consecutive values.
///
/// `erfc` is flat with small arguments (its result changes much less than an
/// ULP between consecutive arguments), so the arguments closest to the inverse
/// of a midpoint have such results (see `values::midpoint_inverse`).
///
/// `erf` with large arguments and `erfc` with large negative arguments are
/// also flat, but all those arguments are tested. Otherwise, the functions
/// are not flat, and the values tested at each binade are dense enough to
/// include results very close to a midpoint.
fn erfc_near_midpoints() -> Vec<f32> {
    let mut xs = Vec::new();
    for x0 in values::binades(-25..=-6, 5000) {
        // erfc(x) = y when erf(x) = 1 - y, which is exact
        let inv = |y: f64| Func::Erf.solve(1.0 - y, x0.into());
        let x = values::midpoint_inverse(x0, fpmath::erfc, inv);
        xs.extend(values::around(x, 1));

        // erfc(-x) = y when erf(x) = y - 1, which is exact
        let inv = |y: f64| Func::Erf.solve(y - 1.0, x0.into());
        let x = values::midpoint_inverse(x0, |x| fpmath::erfc(-x), inv);
        xs.extend(values::around(x, 1));
    }

    // Close arguments can give the same midpoint
    xs.sort_by(f32::total_cmp);
    xs.dedup();
    xs
}

/// Arguments from which `erfcx` is calculated with its asymptotic series.
const ERFCX_TBL_MAX: f32 = 28.0;

/// Arguments from which `erfcx` (of the negated argument) overflows.
const ERFCX_INF_MAX: f32 = 9.5;

fn test_erfcx_with(mut f: impl FnMut(f32)) {
    // Both signs of a value
    let both = |x: f32| [x, -x];

    // Test special values
    for x in values::specials().flat_map(both) {
        f(x);
    }

    // Test a set of subnormal values, whose results are one
    for x in values::subnormals(100_000).flat_map(both) {
        f(x);
    }

    // Test a set of normal values at each binade. The results of the tiny
    // ones are one, the results of the large negative ones overflow, and the
    // results of the largest positive ones are subnormal.
    for x in values::binades(-126..=127, 1000).flat_map(both) {
        f(x);
    }

    // Test more values with small magnitudes
    for x in values::binades(-30..=-2, 200_000).flat_map(both) {
        f(x);
    }

    // Exhaustive test of all the arguments from 0.5 up to where the
    // asymptotic series is used and to the overflow with negative arguments
    for x in values::binades_full(-1..=4) {
        f(x);
    }
    for x in values::interval_full(0.5, ERFCX_INF_MAX) {
        f(-x);
    }

    // Test more values where the polynomial of the asymptotic series is not
    // constant
    for x in values::binades(5..=30, 200_000) {
        f(x);
    }

    // Exhaustive test of all the arguments with subnormal results
    for x in values::interval_full(erfcx_subnormal_min(), f32::MAX) {
        f(x);
    }

    // Test the arguments around the limits of the range
    for x in erfcx_limits() {
        for x in values::around(x, 100_000) {
            f(x);
        }
    }

    // Test the arguments around the boundaries of the argument reduction of
    // the exponential with small arguments (the other ones are all tested)
    for x in reduction_boundaries(0.5).flat_map(both) {
        for x in values::around(x, 2) {
            f(x);
        }
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round
    for x in erfcx_near_midpoints() {
        f(x);
    }
}

/// Returns the smallest argument whose result with `erfcx` is subnormal
/// (`erfcx(x) ~= 1 / (sqrt(π) * x) = 2^-126`).
fn erfcx_subnormal_min() -> f32 {
    (fpmath::scalbn(1.0f64, 126) / std::f64::consts::PI.sqrt()) as f32
}

/// Returns the arguments at the limits of the range of `erfcx`: where it
/// changes the evaluation method (also with negative arguments), where it
/// rounds to one with tiny arguments (positive and negative), where it
/// becomes subnormal and where it overflows.
fn erfcx_limits() -> [f32; 10] {
    let exp2i = |e: i32| fpmath::scalbn(1.0f64, e);
    // erfcx(x) ~= 1 - x * 2 / sqrt(π) with tiny x
    let tiny = |e: i32| (std::f64::consts::PI.sqrt() * exp2i(e - 1)) as f32;
    // 2 * exp(x^2) = 2^128
    let overflow = (127.0 * std::f64::consts::LN_2).sqrt() as f32;
    [
        0.5,
        ERFCX_TBL_MAX,
        erfcx_subnormal_min(),
        tiny(-25),
        tiny(-24),
        -tiny(-25),
        -tiny(-24),
        -0.5,
        -overflow,
        -ERFCX_INF_MAX,
    ]
}

/// Returns arguments whose results with `erfcx` are very close to a midpoint
/// between two consecutive values.
///
/// `erfcx` is flat with small arguments (its result changes much less than an
/// ULP between consecutive arguments), so the arguments closest to the inverse
/// of a midpoint have such results (see `values::midpoint_inverse`).
///
/// Otherwise, it is not flat, and the values tested at each binade (all of
/// them from 0.5) are dense enough to include results very close to a
/// midpoint.
fn erfcx_near_midpoints() -> Vec<f32> {
    // Solves `erfcx(x) = y` in `f64` with Newton's method, starting from
    // `x0`, with erfcx'(x) = 2 * x * erfcx(x) - 2 / sqrt(π)
    let solve = |y: f64, x0: f64| {
        let mut x = x0;
        for _ in 0..100 {
            let fx = fpmath::erfcx(x);
            let dx = (fx - y) / (2.0 * x * fx - std::f64::consts::FRAC_2_SQRT_PI);
            x -= dx;
            if dx.abs() <= x.abs() * f64::EPSILON {
                break;
            }
        }
        x
    };

    let mut xs = Vec::new();
    for x0 in values::binades(-25..=-6, 5000).flat_map(|x| [x, -x]) {
        let x = values::midpoint_inverse(x0, fpmath::erfcx, |y| solve(y, x0.into()));
        xs.extend(values::around(x, 1));
    }

    // Close arguments can give the same midpoint
    xs.sort_by(f32::total_cmp);
    xs.dedup();
    xs
}
