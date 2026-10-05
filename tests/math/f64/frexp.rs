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
            let big = rug::Float::with_val(53, x);
            assert_eq!(e, big.get_exp().unwrap(), "x = {x:?}");
            assert_result_eq!(m, (big >> e).to_f64());
        }
    });
}

fn test_with(mut f: impl FnMut(f64)) {
    // Test special values
    for x in values::specials() {
        f(x);
        f(-x);
    }

    // Test a set of subnormal values, including all the bit lengths (they
    // are normalized differently)
    for x in values::subnormals(100_000) {
        f(x);
        f(-x);
    }

    // Test a set of normal values at each binade
    for x in values::binades(-1022..=1023, 1000) {
        f(x);
        f(-x);
    }
}
