use super::{MIN_MAX_ERROR, check_result, values};

#[test]
fn test_sinh_cosh() {
    let mut max_sin1_error: f32 = 0.0;
    let mut max_sin2_error: f32 = 0.0;
    let mut max_cos1_error: f32 = 0.0;
    let mut max_cos2_error: f32 = 0.0;
    test_with(|x| {
        let (expected_sin, expected_cos) = fpmath::sinh_cosh(f64::from(x));

        let actual_sin1 = fpmath::sinh(x);
        let actual_cos1 = fpmath::cosh(x);
        let (actual_sin2, actual_cos2) = fpmath::sinh_cosh(x);
        assert_result_eq!(fpmath::sinh(-x), -actual_sin1);
        assert_result_eq!(fpmath::cosh(-x), actual_cos1);
        assert_result_eq!(fpmath::sinh_cosh(-x), (-actual_sin2, actual_cos2));

        check_result(x, actual_sin1, expected_sin, &mut max_sin1_error);
        check_result(x, actual_sin2, expected_sin, &mut max_sin2_error);
        check_result(x, actual_cos1, expected_cos, &mut max_cos1_error);
        check_result(x, actual_cos2, expected_cos, &mut max_cos2_error);
    });
    eprintln!("max sinh1 error = {max_sin1_error}");
    eprintln!("max sinh2 error = {max_sin2_error}");
    eprintln!("max cosh1 error = {max_cos1_error}");
    eprintln!("max cosh2 error = {max_cos2_error}");
    assert!(max_sin1_error > MIN_MAX_ERROR);
    assert!(max_sin2_error > MIN_MAX_ERROR);
    assert!(max_cos1_error > MIN_MAX_ERROR);
    assert!(max_cos2_error > MIN_MAX_ERROR);
}

#[test]
fn test_tanh() {
    let mut max_error: f32 = 0.0;
    test_with(|x| {
        let expected = fpmath::tanh(f64::from(x));
        let actual = fpmath::tanh(x);
        assert_result_eq!(fpmath::tanh(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

fn test_with(mut f: impl FnMut(f32)) {
    // Negative arguments are tested by the callers.

    // Test special values
    for x in values::specials() {
        f(x);
    }

    // Test a set of subnormal values, whose results are trivial (`x` or one)
    for x in values::subnormals(10_000) {
        f(x);
    }

    // Test a set of normal values at each binade, most of them have trivial
    // results (`x`, one or infinity)
    for x in values::binades(-126..=127, 1000) {
        f(x);
    }

    // Test more values with small magnitudes, whose results are close to `x`
    // or one
    for x in values::binades(-14..=-1, 200_000) {
        f(x);
    }

    // Test values evenly spread over the whole range of non-trivial results,
    // so all the result binades are tested evenly
    for x in values::interval(0.0, 90.0, 20_000_000) {
        f(x);
    }

    // Test the arguments around the limits of the range
    for x in limits() {
        for x in values::around(x, 100_000) {
            f(x);
        }
    }

    // Test the arguments around the boundaries of the argument reduction
    for x in reduction_boundaries() {
        for x in values::around(x, 2) {
            f(x);
        }
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round
    for x in near_midpoints() {
        f(x);
    }
}

/// Returns the arguments at the limits of the range of non-trivial results:
/// above them, `sinh` and `cosh` overflow (their results are above the
/// midpoint between the largest finite value and infinity), and `tanh` is
/// rounded to one (its result is above the midpoint between one and the
/// previous value).
fn limits() -> [f32; 3] {
    let overflow = fpmath::scalbn(1.0, 128) - fpmath::scalbn(1.0, 103);
    [
        fpmath::asinh(overflow) as f32,
        fpmath::acosh(overflow) as f32,
        fpmath::atanh(1.0 - fpmath::scalbn(1.0, -25)) as f32,
    ]
}

/// Returns the arguments `k * ln(2)` and `(k + 1/2) * ln(2)` (where `k`
/// changes in the argument reduction) up to the overflow limit.
fn reduction_boundaries() -> impl Iterator<Item = f32> {
    (0..=2 * 129).map(|k2| (f64::from(k2) * std::f64::consts::LN_2 / 2.0) as f32)
}

/// Returns arguments whose results are very close to a midpoint between two
/// consecutive values.
///
/// With small arguments, `sinh(x) = x + g(x)` and `tanh(x) = x + g(x)` are not
/// flat, but `g` is (see `values::offset_midpoint`).
///
/// `cosh` is flat close to zero and `tanh` is flat with large arguments
/// (their results change much less than an ULP between consecutive arguments,
/// and they are close to one), so the arguments closest to the inverse of a
/// midpoint have such results (see `values::midpoint_inverse`).
fn near_midpoints() -> Vec<f32> {
    let mut xs = Vec::new();
    for x0 in values::binades(-12..=-2, 20_000) {
        // g(x) = sinh(x) - x, g'(x) = cosh(x) - 1
        let dg = |x: f64| fpmath::cosh(x) - 1.0;
        if let Some(x) = values::offset_midpoint(x0, fpmath::sinh, dg) {
            xs.extend(values::around(x, 1));
        }
        // g(x) = tanh(x) - x, g'(x) = 1 / cosh(x)^2 - 1 = -tanh(x)^2
        let dg = |x: f64| -(fpmath::tanh(x) * fpmath::tanh(x));
        if let Some(x) = values::offset_midpoint(x0, fpmath::tanh, dg) {
            xs.extend(values::around(x, 1));
        }
        xs.extend(values::around(
            values::midpoint_inverse(x0, fpmath::cosh, fpmath::acosh),
            1,
        ));
    }
    for x0 in values::interval(2.5, 9.0, 100_000) {
        xs.extend(values::around(
            values::midpoint_inverse(x0, fpmath::tanh, fpmath::atanh),
            1,
        ));
    }

    // Close arguments can give the same midpoint
    xs.sort_by(f32::total_cmp);
    xs.dedup();
    xs
}
