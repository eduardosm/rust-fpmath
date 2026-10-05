use super::{MIN_MAX_ERROR, check_result, values};

#[test]
fn test_asinh() {
    let mut max_error: f32 = 0.0;
    test_asinh_with(|x| {
        let expected = fpmath::asinh(f64::from(x));
        let actual = fpmath::asinh(x);
        assert_result_eq!(fpmath::asinh(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_acosh() {
    let mut max_error: f32 = 0.0;
    test_acosh_with(|x| {
        let expected = fpmath::acosh(f64::from(x));
        let actual = fpmath::acosh(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_atanh() {
    let mut max_error: f32 = 0.0;
    test_atanh_with(|x| {
        let expected = fpmath::atanh(f64::from(x));
        let actual = fpmath::atanh(x);
        assert_result_eq!(fpmath::atanh(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

fn test_asinh_with(mut f: impl FnMut(f32)) {
    // Negative arguments are tested by the callers.

    // Test special values
    for x in values::specials() {
        f(x);
    }

    // Test a set of subnormal values, whose results are trivial (`x`)
    for x in values::subnormals(10_000) {
        f(x);
    }

    // Test across a set of normal values at each binade
    for x in values::binades(-126..=127, 10_000) {
        f(x);
    }

    // Test more values at the binades with moderate magnitudes
    for x in values::binades(-13..=13, 200_000) {
        f(x);
    }

    // Test the arguments around the limit below which the implementation
    // returns `x` (2^-20)
    for x in values::around(fpmath::scalbn(1.0, -20), 100_000) {
        f(x);
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round
    for x in asinh_near_midpoints() {
        f(x);
    }
}

fn test_acosh_with(mut f: impl FnMut(f32)) {
    // Test special values
    for x in values::specials() {
        f(x);
        f(-x);
    }

    // Test across a set of normal values at each binade above one, and the
    // mantissa patterns below one and of negative values (whose results are
    // NaN)
    for x in values::binades(0..=127, 20_000) {
        f(x);
    }
    for x in values::binades(-126..=-1, 0).chain(values::subnormals(0)) {
        f(x);
        f(-x);
    }
    for x in values::binades(0..=127, 0) {
        f(-x);
    }

    // Exhaustive test of all the values in [1, 2], whose results close to one
    // are close to zero, and some values below one (whose results are NaN)
    for x in values::interval_full(1.0, 2.0).chain(values::around(1.0, 100)) {
        f(x);
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round. With large arguments, the
    // function is flat (its result changes much less than an ULP between
    // consecutive arguments), so the arguments closest to the inverse of a
    // midpoint have such results (see `values::midpoint_inverse`).
    for x0 in values::binades(16..=127, 1000) {
        let x = values::midpoint_inverse(x0, fpmath::acosh, fpmath::cosh);
        for x in values::around(x, 1) {
            f(x);
        }
    }
}

fn test_atanh_with(mut f: impl FnMut(f32)) {
    // Negative arguments are tested by the callers.

    // Test special values
    for x in values::specials() {
        f(x);
    }

    // Test a set of subnormal values, whose results are trivial (`x`)
    for x in values::subnormals(10_000) {
        f(x);
    }

    // Test across a set of normal values at each binade below one, and the
    // mantissa patterns at the binades above one (whose results are NaN)
    for x in values::binades(-126..=-1, 20_000).chain(values::binades(0..=127, 0)) {
        f(x);
    }

    // Test more values at the binades where the results are not just `x`
    for x in values::binades(-13..=-2, 200_000) {
        f(x);
    }

    // Exhaustive test of all the values in [1/2, 1], whose results close to
    // one are large, and some values above one (whose results are NaN)
    for x in values::interval_full(0.5, 1.0).chain(values::around(1.0, 100)) {
        f(x);
    }

    // Test the arguments around the limit below which the implementation
    // returns `x` (2^-20)
    for x in values::around(fpmath::scalbn(1.0, -20), 100_000) {
        f(x);
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round. With small arguments,
    // `atanh(x) = x + g(x)` is not flat, but `g` is (see
    // `values::offset_midpoint`).
    // g(x) = atanh(x) - x, g'(x) = 1 / (1 - x^2) - 1 = x^2 / (1 - x^2)
    let dg = |x: f64| x * x / (1.0 - x * x);
    for x0 in values::binades(-12..=-2, 20_000) {
        if let Some(x) = values::offset_midpoint(x0, fpmath::atanh, dg) {
            for x in values::around(x, 1) {
                f(x);
            }
        }
    }
}

/// Returns arguments whose results with `asinh` are very close to a midpoint
/// between two consecutive values.
///
/// With small arguments, `asinh(x) = x + g(x)` is not flat, but `g` is (see
/// `values::offset_midpoint`). With large arguments, `asinh` is flat (its
/// result changes much less than an ULP between consecutive arguments), so
/// the arguments closest to the inverse of a midpoint have such results (see
/// `values::midpoint_inverse`).
fn asinh_near_midpoints() -> Vec<f32> {
    let mut xs = Vec::new();

    // g(x) = asinh(x) - x, g'(x) = 1 / sqrt(1 + x^2) - 1
    let dg = |x: f64| 1.0 / (1.0 + x * x).sqrt() - 1.0;
    for x0 in values::binades(-12..=-2, 20_000) {
        if let Some(x) = values::offset_midpoint(x0, fpmath::asinh, dg) {
            xs.extend(values::around(x, 1));
        }
    }

    for x0 in values::binades(16..=127, 1000) {
        xs.extend(values::around(
            values::midpoint_inverse(x0, fpmath::asinh, fpmath::sinh),
            1,
        ));
    }

    // Close arguments can give the same midpoint
    xs.sort_by(f32::total_cmp);
    xs.dedup();
    xs
}
