use super::values;

#[test]
fn test_scalbn() {
    test_with(|x, n| {
        let actual = fpmath::scalbn(x, n);
        assert_result_eq!(fpmath::scalbn(-x, n), -actual);

        // `x * 2^n` is exact in `f64` when it is in its range, and otherwise
        // the result is zero or infinity in both types, so it is rounded only
        // once to `f32`
        let expected = fpmath::scalbn(f64::from(x), n) as f32;
        assert_result_eq!(actual, expected);
    });
}

fn test_with(mut f: impl FnMut(f32, i32)) {
    // Negative arguments are tested by the caller.

    // Test special values with special exponents (including the ones that
    // move values between the limits of the range, and the extreme ones)
    let ns = [
        0,
        1,
        2,
        126,
        127,
        128,
        149,
        150,
        252,
        253,
        254,
        275,
        276,
        277,
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
    let exps = [-126, -125, -1, 0, 1, 126, 127];
    for x in values::binades(exps, 1000).chain(values::subnormals(1000)) {
        // Exponent of `x`
        let e = rug::Float::with_val(24, x).get_exp().unwrap() - 1;
        for t in -155..=132 {
            f(x, t - e);
        }
    }

    // Test the rounding of subnormal results: a dense set of mantissas
    // (including the patterns that give ties and carries) with results at
    // every exponent from below half the smallest subnormal to the smallest
    // normal
    for x in values::binades([0], 1_000_000) {
        for n in -152..=-125 {
            f(x, n);
        }
    }

    // Test a set of values at each binade with exponents evenly spread from
    // -300 to 300, which cover all the results
    for x in values::binades(-126..=127, 100) {
        for u in crate::utils::spread(32, 20, x.to_bits().into()) {
            f(x, (u % 601) as i32 - 300);
        }
    }
}
