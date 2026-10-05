use super::values;

#[test]
fn test_frexp() {
    test_with(|x| {
        let (m, e) = fpmath::frexp(x);
        if x == 0.0 || !x.is_finite() {
            assert_result_eq!(m, x);
            assert_eq!(e, 0);
        } else {
            // `x = m * 2^e` with `0.5 <= |m| < 1`, like the exponent of `rug`
            let big = rug::Float::with_val(24, x);
            assert_eq!(e, big.get_exp().unwrap(), "x = {x:?}");
            assert_result_eq!(m, (big >> e).to_f32());
        }
    });
}

fn test_with(mut f: impl FnMut(f32)) {
    // Test special values
    for x in values::specials() {
        f(x);
        f(-x);
    }

    // Exhaustive test of all subnormal numbers (they are normalized
    // differently depending on their bit length)
    for x in values::subnormals_full() {
        f(x);
        f(-x);
    }

    // Test across a wide range of normal numbers
    for x in values::binades(-126..=127, 10_000) {
        f(x);
        f(-x);
    }
}
