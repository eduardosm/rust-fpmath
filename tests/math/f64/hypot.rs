use super::{MIN_MAX_ERROR, RUG_PREC, check_result, exponent, values};

#[test]
fn test_hypot() {
    let mut max_error: f64 = 0.0;
    test_with(|x, y| {
        let expected = rug::Float::with_val(RUG_PREC, x).hypot(&rug::Float::with_val(RUG_PREC, y));
        let actual = fpmath::hypot(x, y);
        assert_result_eq!(fpmath::hypot(-x, y), actual);
        assert_result_eq!(fpmath::hypot(x, -y), actual);
        assert_result_eq!(fpmath::hypot(-x, -y), actual);
        assert_result_eq!(fpmath::hypot(y, x), actual);
        assert_result_eq!(fpmath::hypot(-y, x), actual);
        assert_result_eq!(fpmath::hypot(y, -x), actual);
        assert_result_eq!(fpmath::hypot(-y, -x), actual);

        check_result((x, y), actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

fn test_with(mut f: impl FnMut(f64, f64)) {
    // Toggling the signs of `x` and `y` is done by the caller.

    // Test special values
    for x in values::specials() {
        for y in values::specials() {
            f(x, y);
        }
    }

    // Test with zero
    for x in values::binades(-1022..=1023, 10).chain(values::subnormals(100_000)) {
        f(x, 0.0);
    }

    // Test the exponent differences where both arguments affect the result
    // (the smaller one stops affecting it when it is less than 2^-27 times
    // the larger one), at every scale, including subnormals
    for (x, y) in values::binade_pairs(values::exp_pairs(0..=30, 1), 2) {
        f(x, y);
    }

    // Test larger exponent differences at some scales
    for (x, y) in values::binade_pairs(values::exp_pairs((31..=2097).step_by(7), 29), 1) {
        f(x, y);
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round
    for (x, y) in near_midpoints() {
        f(x, y);
    }

    // Test exact results
    for (x, y) in pythagorean_triples() {
        f(x, y);
    }
}

/// Returns pairs `(a, b)` such that `hypot(a, b)` is very close to a midpoint
/// between two consecutive values.
///
/// With `0 < b < a`, `hypot(a, b) = a + (k + 1/2) * ulp(a)` when
/// `b = sqrt((a + (k + 1/2) * ulp(a))^2 - a^2) ~= a * sqrt((2 * k + 1) * 2^-52)`.
/// Rounding `b` moves the result by about `(b / a)^2 / 2` ULPs, which is less
/// than 2^-31 ULPs with `k < 2^20`, and the neighbors of `b` move it to both
/// sides of the midpoint. `a` is large enough to keep `b` normal, since a
/// subnormal `b` would be too imprecise.
fn near_midpoints() -> impl Iterator<Item = (f64, f64)> {
    values::binades(-995..=1023, 0).flat_map(|a| {
        // Number of ULPs from `a` to the end of its binade, the midpoint
        // must be before it
        let ulps_left = (1 << 53) - ((a.to_bits() & ((1 << 52) - 1)) | (1 << 52));
        [0, 1, 2]
            .into_iter()
            .chain(crate::utils::spread(20, 4, a.to_bits()))
            .map(move |k| k % ulps_left)
            .flat_map(move |k| {
                const PREC: u32 = 256;
                let a_big = rug::Float::with_val(PREC, a);
                let half_ulps = rug::Float::with_val(PREC, 2 * k + 1) << (exponent(a) - 53);
                let h = rug::Float::with_val(PREC, &a_big + &half_ulps);
                let b2 = rug::Float::with_val(PREC, &h - &a_big) * (h + a_big);
                let b = b2.sqrt().to_f64();
                values::around(b, 1).map(move |b| (a, b))
            })
    })
}

/// Returns pairs `(a, b)` such that `hypot(a, b)` is exact, scaled by powers of
/// two so that the result takes every exponent.
fn pythagorean_triples() -> impl Iterator<Item = (f64, f64)> {
    // Small triples, and large triples whose results have up to 53
    // significant bits
    let small = (2..=32).flat_map(|m| (1..m).map(move |n| (m, n)));
    let large = crate::utils::spread2(25, 64, 0).map(|(m, n)| ((1 << 24) + m, 1 + (n >> 1)));

    small.chain(large).flat_map(|(m, n): (u64, u64)| {
        // Euclid's formula: (m^2 - n^2)^2 + (2 * m * n)^2 = (m^2 + n^2)^2
        let a = (m * m - n * n) as f64;
        let b = (2 * m * n) as f64;
        let c = (m * m + n * n) as f64;
        (-1074..=1023).filter_map(move |e| {
            let s = e - exponent(c);
            let scale = |v: f64| {
                let big = rug::Float::with_val(64, v) << s;
                let r = big.to_f64();
                (big == r).then_some(r)
            };
            // Skip when the arguments or the result are not exact
            scale(c)?;
            Some((scale(a)?, scale(b)?))
        })
    })
}
