use super::{MIN_MAX_ERROR, check_result, exponent, values};

#[test]
fn test_hypot() {
    let mut max_error: f32 = 0.0;
    test_with(|x, y| {
        let expected = fpmath::hypot(f64::from(x), f64::from(y));
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

fn test_with(mut f: impl FnMut(f32, f32)) {
    // Toggling the signs of `x` and `y` is done by the caller.

    // Test special values
    for x in values::specials() {
        for y in values::specials() {
            f(x, y);
        }
    }

    // Test with zero
    for x in values::binades(-126..=127, 10).chain(values::subnormals_full()) {
        f(x, 0.0);
    }

    // Test the exponent differences where both arguments affect the result
    // (the smaller one stops affecting it when it is less than 2^-12 times
    // the larger one), at every scale, including subnormals
    for (x, y) in values::binade_pairs(values::exp_pairs(0..=30, 1), 200) {
        f(x, y);
    }

    // Test larger exponent differences at every scale
    for (x, y) in values::binade_pairs(values::exp_pairs(31..=276, 1), 16) {
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
/// `b = sqrt((a + (k + 1/2) * ulp(a))^2 - a^2) ~= a * sqrt((2 * k + 1) * 2^-23)`.
/// Rounding `b` moves the result by about `(b / a)^2 / 2` ULPs, which is less
/// than 2^-12 ULPs with `k < 2^10`, and the neighbors of `b` move it to both
/// sides of the midpoint. `a` is large enough to keep `b` normal, since a
/// subnormal `b` would be too imprecise.
fn near_midpoints() -> impl Iterator<Item = (f32, f32)> {
    values::binades(-114..=127, 20).flat_map(|a| {
        // Number of ULPs from `a` to the end of its binade, the midpoint
        // must be before it
        let ulps_left = (1 << 24) - ((a.to_bits() & ((1 << 23) - 1)) | (1 << 23));
        [0, 1, 2]
            .into_iter()
            .chain(crate::utils::spread(10, 8, a.to_bits().into()).map(|k| k as u32))
            .map(move |k| k % ulps_left)
            .flat_map(move |k| {
                // All the operations are exact in f64, except the square root
                let a64 = f64::from(a);
                let half_ulps = fpmath::scalbn(f64::from(2 * k + 1), exponent(a) - 24);
                let h = a64 + half_ulps;
                let b = ((h - a64) * (h + a64)).sqrt() as f32;
                values::around(b, 1).map(move |b| (a, b))
            })
    })
}

/// Returns pairs `(a, b)` such that `hypot(a, b)` is exact, scaled by powers of
/// two so that the result takes every exponent.
fn pythagorean_triples() -> impl Iterator<Item = (f32, f32)> {
    // Small triples, and large triples whose results have up to 24
    // significant bits
    let small = (2..=32).flat_map(|m| (1..m).map(move |n| (m, n)));
    let large =
        crate::utils::spread2(11, 256, 0).map(|(m, n)| ((1 << 10) + m as u32, 1 + (n >> 1) as u32));

    small.chain(large).flat_map(|(m, n): (u32, u32)| {
        // Euclid's formula: (m^2 - n^2)^2 + (2 * m * n)^2 = (m^2 + n^2)^2
        let a = (m * m - n * n) as f32;
        let b = (2 * m * n) as f32;
        let c = (m * m + n * n) as f32;
        (-149..=127).filter_map(move |e| {
            let scale = |v: f32| {
                let big = fpmath::scalbn(f64::from(v), e - exponent(c));
                let r = big as f32;
                (f64::from(r) == big).then_some(r)
            };
            // Skip when the arguments or the result are not exact
            scale(c)?;
            Some((scale(a)?, scale(b)?))
        })
    })
}
