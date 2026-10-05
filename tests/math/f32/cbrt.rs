use super::{MIN_MAX_ERROR, check_result, values};

#[test]
fn test_cbrt() {
    let mut max_error: f32 = 0.0;
    test_with(|x| {
        let expected = fpmath::cbrt(f64::from(x));
        let actual = fpmath::cbrt(x);
        assert_result_eq!(fpmath::cbrt(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

fn test_with(mut f: impl FnMut(f32)) {
    // Negative arguments are tested by the caller.

    // Test special values
    for x in values::specials() {
        f(x);
    }

    // Exhaustive test of all subnormal numbers
    for x in values::subnormals_full() {
        f(x);
    }

    // Exhaustive test of all mantissas at three consecutive binades. The
    // implementation only depends on the mantissa and on the exponent modulo 3
    // (the result is scaled exactly), so this covers all the arguments.
    for x in values::binades_full(0..=2) {
        f(x);
    }

    // Test across a wide range of normal numbers
    for x in values::binades(-126..=127, 500) {
        f(x);
    }

    // Test integers
    for arg in 1..=500_000 {
        f(arg as f32);
    }
}
