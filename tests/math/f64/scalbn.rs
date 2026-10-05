use super::values;

#[test]
fn test_scalbn() {
    test_with(|x, n| {
        let actual = fpmath::scalbn(x, n);
        assert_result_eq!(fpmath::scalbn(-x, n), -actual);

        // `x * 2^n` is exact with `rug`, and then it is rounded
        let expected = (rug::Float::with_val(53, x) << n).to_f64();
        assert_result_eq!(actual, expected);
    });
}

fn test_with(mut f: impl FnMut(f64, i32)) {
    // Negative arguments are tested by the caller.

    // Test special values with special exponents (including the ones that
    // move values between the limits of the range, and the extreme ones)
    let ns = [
        0,
        1,
        2,
        1022,
        1023,
        1024,
        1074,
        1075,
        2045,
        2046,
        2047,
        2097,
        2098,
        2099,
        i32::MAX - 1,
        i32::MAX,
    ];
    for x in values::specials().chain([1.0]) {
        for n in ns {
            f(x, n);
            f(x, -n);
        }
        f(x, i32::MIN);
    }

    // Test the results at every exponent, from beyond the overflow limit to
    // below half the smallest subnormal, for a set of values at the binades at
    // the limits and in the middle of the range (and subnormals). The exact
    // results only depend on the mantissas and on their exponents, which
    // determine how subnormal results are rounded.
    let exps = [-1022, -1021, -1, 0, 1, 1022, 1023];
    for x in values::binades(exps, 100).chain(values::subnormals(100)) {
        // Exponent of `x`
        let e = rug::Float::with_val(53, x).get_exp().unwrap() - 1;
        for t in -1080..=1030 {
            f(x, t - e);
        }
    }

    // Test a set of values at each binade with exponents evenly spread from
    // -2200 to 2200, which cover all the results
    for x in values::binades(-1022..=1023, 10) {
        for u in crate::utils::spread(32, 10, x.to_bits()) {
            f(x, (u % 4401) as i32 - 2200);
        }
    }
}
