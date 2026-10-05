use super::{MIN_MAX_ERROR, RUG_PREC, check_result, values};

#[test]
fn test_cbrt() {
    let mut max_error: f64 = 0.0;
    test_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).cbrt();
        let actual = fpmath::cbrt(x);
        assert_result_eq!(fpmath::cbrt(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

fn test_with(mut f: impl FnMut(f64)) {
    // Negative arguments are tested by the caller.

    // Test special values
    for x in values::specials() {
        f(x);
    }

    // Test a set of subnormal values
    for x in values::subnormals(1_000_000) {
        f(x);
    }

    // Test across a set of normal values at each binade
    for x in values::binades(-1022..=1023, 1000) {
        f(x);
    }

    // Test more values at three consecutive binades. The implementation only
    // depends on the mantissa and on the exponent modulo 3 (the result is
    // scaled exactly).
    for x in values::binades(0..=2, 2_000_000) {
        f(x);
    }

    // Test exact results (and their neighbors) at every binade
    for x in exact_cubes() {
        for x in values::around(x, 1) {
            f(x);
        }
    }

    // Test integers
    for arg in 1..=500_000 {
        f(arg as f64);
    }
}

/// Returns arguments with exact cube roots: cubes of integers in
/// `[2^(50/3), 2^(53/3))` (which are exact and cover three binades) scaled by
/// powers of eight.
fn exact_cubes() -> impl Iterator<Item = f64> {
    let lo = 104_032; // ceil(2^(50/3))
    let hi = 208_063; // floor(2^(53/3))
    crate::utils::spread(32, 200, 0)
        .map(move |u| lo + u % (hi - lo + 1))
        .flat_map(|r| {
            let x = (r * r * r) as f64;
            (-358..=324).map(move |k| fpmath::scalbn(x, 3 * k))
        })
        .filter(|x| x.is_normal())
}
