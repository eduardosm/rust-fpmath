use super::values;

#[test]
fn test_round() {
    test_with(|arg| {
        let expected = fpmath::round(f64::from(arg)) as f32;
        let actual = fpmath::round(arg);
        assert_result_eq!(expected, actual);
    });
}

#[test]
fn test_floor() {
    test_with(|arg| {
        let expected = fpmath::floor(f64::from(arg)) as f32;
        let actual = fpmath::floor(arg);
        assert_result_eq!(expected, actual);
    });
}

#[test]
fn test_ceil() {
    test_with(|arg| {
        let expected = fpmath::ceil(f64::from(arg)) as f32;
        let actual = fpmath::ceil(arg);
        assert_result_eq!(expected, actual);
    });
}

#[test]
fn test_trunc() {
    test_with(|arg| {
        let expected = fpmath::trunc(f64::from(arg)) as f32;
        let actual = fpmath::trunc(arg);
        assert_result_eq!(expected, actual);
    });
}

fn test_with(f: fn(f32)) {
    // Test special values
    for x in values::specials() {
        f(x);
        f(-x);
    }

    // Exhaustive test of all subnormal numbers, whose results are zero or
    // one (with the sign of the argument)
    for x in values::subnormals_full() {
        f(x);
        f(-x);
    }

    // Test across a wide range of normal numbers
    for x in values::binades(-126..=127, 10_000) {
        f(x);
        f(-x);
    }

    // Exhaustive test of all the values from 1/2 to 2^23, whose results
    // depend on their fractional parts (smaller values are rounded to zero or
    // one, and larger values are integers)
    for x in values::binades_full(-1..=22) {
        f(x);
        f(-x);
    }
}
