use super::values;

#[test]
fn test_sqrt() {
    test_with(|x| {
        assert!(x.partial_cmp(&0.0) != Some(std::cmp::Ordering::Less));
        let actual = fpmath::sqrt(x);
        let expected = rug::Float::with_val(24, x).sqrt().to_f32();
        if x != 0.0 {
            assert!(fpmath::sqrt(-x).is_nan());
        }
        assert_result_eq!(actual, expected);
    });
}

fn test_with(mut f: impl FnMut(f32)) {
    // Test special values
    for x in values::specials().chain([-0.0]) {
        f(x);
    }

    // Exhaustive test of all subnormal numbers
    for x in values::subnormals_full() {
        f(x);
    }

    // Exhaustive test of all mantissas at two consecutive binades. The
    // implementation only depends on the mantissa and on the parity of the
    // exponent (the result is scaled exactly), so this covers all the
    // arguments.
    for x in values::binades_full(0..=1) {
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
