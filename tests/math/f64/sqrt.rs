use super::values;

#[test]
fn test_sqrt() {
    test_with(|x| {
        let actual = fpmath::sqrt(x);
        let expected = rug::Float::with_val(53, x).sqrt().to_f64();
        if x != 0.0 {
            assert!(fpmath::sqrt(-x).is_nan());
        }
        assert_result_eq!(actual, expected);
    });
}

fn test_with(mut f: impl FnMut(f64)) {
    // Test special values
    for x in values::specials().chain([-0.0]) {
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

    // Test more values at two consecutive binades. The implementation only
    // depends on the mantissa and on the parity of the exponent (the result
    // is scaled exactly).
    for x in values::binades(0..=1, 5_000_000) {
        f(x);
    }

    // Test exact results (and their neighbors) at every binade
    for x in exact_squares() {
        for x in values::around(x, 1) {
            f(x);
        }
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round, at every binade
    for x in near_midpoints() {
        f(x);
    }

    // Test integers
    for arg in 1..=500_000 {
        f(arg as f64);
    }
}

/// Returns arguments with exact square roots: squares of integers in
/// `[2^25.5, 2^26.5)` (which are exact and cover two binades) scaled by powers
/// of four.
fn exact_squares() -> impl Iterator<Item = f64> {
    let lo = 47_453_133; // ceil(2^25.5)
    let hi = 94_906_265; // floor(2^26.5)
    crate::utils::spread(32, 200, 0)
        .map(move |u| lo + u % (hi - lo + 1))
        .flat_map(|r| {
            let x = (r * r) as f64;
            (-537..=486).map(move |k| fpmath::scalbn(x, 2 * k))
        })
        .filter(|x| x.is_normal())
}

/// Returns arguments whose square roots are very close to a midpoint between
/// two consecutive values, at every binade.
///
/// With `x = M * 2^(e - 52)` (`M` a 53-bit integer) and `p` the parity of `e`,
/// `sqrt(x) = sqrt(n) * 2^((e - p) / 2 - 52)` with `n = M * 2^(52 + p)`, where
/// `sqrt(n)` is in `[2^52, 2^53)` and has to be rounded to an integer. When
/// `4 * n = a^2 - d` with an odd `a` and a small `d`,
/// `sqrt(n) = a / 2 - d / (4 * a) + ...`, which is very close to the midpoint
/// `a / 2` (`|d| / 2^55` ULPs at most). Such `a` are the square roots of `d`
/// modulo `2^(54 + p)`, which exist when `d = 1 (mod 8)`.
fn near_midpoints() -> Vec<f64> {
    let mut xs = Vec::new();
    for d in (-50..=50).map(|j| 1 + 8 * j) {
        for p in [0, 1] {
            let k = 54 + p;
            let roots = sqrt_mod_pow2(d, k);
            for a in roots.filter(|a| (1 << 53..1 << 54).contains(a)) {
                // M = (a^2 - d) / 2^k, exact
                let m = (a * a).wrapping_sub(d as u128) >> k;
                if (1 << 52..1 << 53).contains(&m) {
                    let m = m as f64;
                    let e_lo = -1022 + (p as i32 + 1022) % 2;
                    xs.extend((e_lo..=1023).step_by(2).map(|e| fpmath::scalbn(m, e - 52)));
                }
            }
        }
    }
    xs
}

/// Returns the four odd square roots of `d` modulo `2^k`, with
/// `d = 1 (mod 8)`.
fn sqrt_mod_pow2(d: i128, k: u32) -> impl Iterator<Item = u128> {
    assert_eq!(d.rem_euclid(8), 1);
    let mask = |i: u32| (1u128 << i) - 1;
    let d = d as u128 & mask(k);

    // Hensel lifting: if `a^2 = d (mod 2^i)` (with i >= 3), then
    // `a^2 = d (mod 2^(i + 1))` or `(a + 2^(i - 1))^2 = d (mod 2^(i + 1))`
    let mut a = 1u128;
    for i in 3..k {
        if ((a * a) ^ d) & mask(i + 1) != 0 {
            a += 1 << (i - 1);
        }
    }

    // The other roots are `-a` and `±a + 2^(k - 1)`
    let half = 1 << (k - 1);
    [
        a,
        a.wrapping_neg(),
        a + half,
        a.wrapping_neg().wrapping_add(half),
    ]
    .into_iter()
    .map(move |a| a & mask(k))
}
