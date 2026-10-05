use super::{MIN_MAX_ERROR, RUG_PREC, check_result, purify, values};

#[test]
fn test_asinh() {
    let mut max_error: f64 = 0.0;
    test_asinh_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).asinh();
        let actual = fpmath::asinh(x);
        assert_result_eq!(fpmath::asinh(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_acosh() {
    let mut max_error: f64 = 0.0;
    test_acosh_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).acosh();
        let actual = fpmath::acosh(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_atanh() {
    let mut max_error: f64 = 0.0;
    test_atanh_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).atanh();
        let actual = fpmath::atanh(x);
        assert_result_eq!(fpmath::atanh(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

/// Precision used to calculate the test values.
const EXT_PREC: u32 = 256;

fn asinh(x: &rug::Float) -> rug::Float {
    rug::Float::with_val(EXT_PREC, x).asinh()
}

fn acosh(x: &rug::Float) -> rug::Float {
    rug::Float::with_val(EXT_PREC, x).acosh()
}

fn atanh(x: &rug::Float) -> rug::Float {
    rug::Float::with_val(EXT_PREC, x).atanh()
}

fn test_asinh_with(mut f: impl FnMut(f64)) {
    // Negative arguments are tested by the callers.

    // Test special values
    for x in values::specials() {
        f(x);
    }

    // Test a set of subnormal values, whose results are trivial (`x`)
    for x in values::subnormals(10_000) {
        f(x);
    }

    // Test a set of normal values at each binade
    for x in values::binades(-1022..=1023, 500) {
        f(x);
    }

    // Test more values at the binades with moderate magnitudes
    for x in values::binades(-30..=30, 15_000) {
        f(x);
    }

    // Test the arguments around the limits where the implementation changes
    // the evaluation method (2^-6 and 2^501)
    for x in [fpmath::scalbn(1.0, -6), fpmath::scalbn(1.0, 501)] {
        for x in values::around(x, 10_000) {
            f(x);
        }
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round
    for x in asinh_near_midpoints() {
        f(x);
    }
}

fn test_acosh_with(mut f: impl FnMut(f64)) {
    // Test special values
    for x in values::specials() {
        f(x);
        f(-x);
    }

    // Test a set of normal values at each binade above one, and the mantissa
    // patterns below one and of negative values (whose results are NaN)
    for x in values::binades(0..=1023, 2000) {
        f(x);
    }
    for x in values::binades(-1022..=-1, 0).chain(values::subnormals(0)) {
        f(x);
        f(-x);
    }
    for x in values::binades(0..=1023, 0) {
        f(-x);
    }

    // Test all the values within 2^16 ULPs of one, whose results are close to
    // zero (or NaN below one)
    for x in values::around(1.0, 1 << 16) {
        f(x);
    }

    // Test values close to one at larger distances (`1 + v` with `v` at each
    // binade below one)
    for v in values::binades(-36..=-1, 10_000) {
        f(purify(1.0 + v));
    }

    // Test the arguments around the limit where the implementation changes the
    // evaluation method (2^501)
    for x in values::around(fpmath::scalbn(1.0, 501), 10_000) {
        f(x);
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round. With large arguments, the
    // function is flat (its result changes much less than an ULP between
    // consecutive arguments), so the arguments closest to the inverse of a
    // midpoint have such results (see `values::midpoint_inverse`).
    let cosh = |y: &rug::Float| rug::Float::with_val(EXT_PREC, y).cosh();
    for x0 in values::binades((64..=1023).step_by(2), 20) {
        for x in values::around(values::midpoint_inverse(x0, acosh, cosh), 1) {
            f(x);
        }
    }
}

fn test_atanh_with(mut f: impl FnMut(f64)) {
    // Negative arguments are tested by the callers.

    // Test special values
    for x in values::specials() {
        f(x);
    }

    // Test a set of subnormal values, whose results are trivial (`x`)
    for x in values::subnormals(10_000) {
        f(x);
    }

    // Test a set of normal values at each binade below one, and the mantissa
    // patterns at the binades above one (whose results are NaN)
    for x in values::binades(-1022..=-1, 1000).chain(values::binades(0..=1023, 0)) {
        f(x);
    }

    // Test more values at the binades where the results are not just `x`
    for x in values::binades(-30..=-1, 30_000) {
        f(x);
    }

    // Test all the values within 2^16 ULPs of one, whose results are large (or
    // NaN above one)
    for x in values::around(1.0, 1 << 16) {
        f(x);
    }

    // Test values close to one at larger distances (`1 - v` with `v` at each
    // binade below one half)
    for v in values::binades(-37..=-2, 10_000) {
        f(purify(1.0 - v));
    }

    // Test the arguments around the limit where the implementation changes the
    // evaluation method (2^-13)
    for x in values::around(fpmath::scalbn(1.0, -13), 10_000) {
        f(x);
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round. With small arguments,
    // `atanh(x) = x + g(x)` is not flat, but `g` is (see
    // `values::offset_midpoint`).
    // g(x) = atanh(x) - x, g'(x) = 1 / (1 - x^2) - 1 = x^2 / (1 - x^2)
    let dg = |x: &rug::Float| -> rug::Float {
        let x2 = rug::Float::with_val(EXT_PREC, x.square_ref());
        rug::Float::with_val(EXT_PREC, &x2 / (1 - x2.clone()))
    };
    for x0 in values::binades(-26..=-2, 1000) {
        if let Some(x) = values::offset_midpoint(x0, atanh, dg) {
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
fn asinh_near_midpoints() -> Vec<f64> {
    let mut xs = Vec::new();

    // g(x) = asinh(x) - x, g'(x) = 1 / sqrt(1 + x^2) - 1
    let dg = |x: &rug::Float| -> rug::Float {
        let r = rug::Float::with_val(EXT_PREC, 1 + x.clone().square());
        r.sqrt().recip() - 1
    };
    for x0 in values::binades(-26..=-2, 1000) {
        if let Some(x) = values::offset_midpoint(x0, asinh, dg) {
            xs.extend(values::around(x, 1));
        }
    }

    let sinh = |y: &rug::Float| rug::Float::with_val(EXT_PREC, y).sinh();
    for x0 in values::binades((64..=1023).step_by(2), 20) {
        xs.extend(values::around(values::midpoint_inverse(x0, asinh, sinh), 1));
    }

    // Close arguments can give the same midpoint
    xs.sort_by(f64::total_cmp);
    xs.dedup();
    xs
}
