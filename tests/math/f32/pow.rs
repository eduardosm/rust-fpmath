use super::{MIN_MAX_ERROR, check_result, purify, values};

#[test]
fn test_pow() {
    let mut max_error: f32 = 0.0;
    test_pow_with(|x, y| {
        let expected = fpmath::pow(f64::from(x), f64::from(y));
        let actual = fpmath::pow(x, y);
        check_neg_base(x, y);

        check_result((x, y), actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

/// Checks `pow(-|x|, y)` against `pow(|x|, y)`: they are equal when `y` is an
/// even integer (or infinite) and opposite when it is an odd integer. When `y`
/// is not an integer, `pow(-|x|, y)` is NaN if `x` is finite and non-zero, and
/// equal to `pow(|x|, y)` otherwise.
fn check_neg_base(x: f32, y: f32) {
    let pos = fpmath::pow(x.abs(), y);
    let neg = fpmath::pow(-x.abs(), y);
    if y.is_infinite() || y % 2.0 == 0.0 {
        assert_result_eq!(neg, pos);
    } else if y.trunc() == y {
        assert_result_eq!(neg, -pos);
    } else if x.is_finite() && x != 0.0 {
        assert!(neg.is_nan(), "pow({:?}, {y:?}) = {neg:?}", -x.abs());
    } else {
        assert_result_eq!(neg, pos);
    }
}

#[test]
fn test_powi() {
    let mut max_error: f32 = 0.0;
    test_powi_with(|x, y| {
        let expected = fpmath::powi(f64::from(x), y);
        let actual = fpmath::powi(x, y);
        if f64::from(y as f32) == f64::from(y) {
            assert_result_eq!(fpmath::pow(x, y as f32), actual);
        }
        if y % 2 == 0 {
            assert_result_eq!(fpmath::powi(-x, y), actual);
        } else {
            assert_result_eq!(fpmath::powi(-x, y), -actual);
        }

        check_result((x, y), actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

fn test_pow_with(mut f: impl FnMut(f32, f32)) {
    // Negative bases are tested by `check_neg_base`.

    // Test special values
    for x in values::specials().chain([1.0]) {
        for y in values::specials().chain(special_exponents()) {
            f(x, y);
            f(x, -y);
        }
    }

    // Test halves and small integers as exponents with bases at every binade
    for x in values::binades(-126..=127, 100).chain(values::subnormals(1000)) {
        for y in [0.5, 1.0, 1.5, 2.0, 3.0] {
            f(x, y);
            f(x, -y);
        }
    }

    // Test the largest exponents with bases close to one
    for x in values::binades([-1, 0], 1000) {
        for y in special_exponents() {
            f(x, y);
            f(x, -y);
        }
    }

    // Test results spread across the whole range, with bases at every binade,
    // including subnormals and bases close to one
    for x in test_bases() {
        for y in result_targets(x) {
            f(x, y);
        }
    }

    // Test exact results and ties
    for (x, y) in exact_results() {
        f(x, y);
    }
}

fn test_powi_with(mut f: impl FnMut(f32, i32)) {
    // Negative bases are tested by the caller.

    // Test special values
    for x in values::specials().chain([1.0]) {
        for y in [0, 1, 2, 3, i32::MAX - 1, i32::MAX] {
            f(x, y);
            f(x, -y);
        }
        f(x, i32::MIN);
    }

    // Test small exponents with bases at every binade
    for x in values::binades(-126..=127, 100).chain(values::subnormals(1000)) {
        for y in [1, 2, 3] {
            f(x, y);
            f(x, -y);
        }
    }

    // Test results spread across the whole range, with bases at every binade,
    // including subnormals and bases close to one (where the exponents can
    // saturate)
    for x in test_bases() {
        for y in result_targets(x) {
            if y.trunc() == y {
                f(x, y as i32);
            }
        }
    }

    // Test exact results and ties
    for (x, y) in exact_results() {
        if y.trunc() == y {
            f(x, y as i32);
        }
    }
}

/// Returns positive exponents that need special handling: halves and small
/// integers, and the largest integers and non-integers, where the parity
/// matters for negative bases.
fn special_exponents() -> impl Iterator<Item = f32> {
    const P23: f32 = (1u32 << 23) as f32;
    const P24: f32 = (1u32 << 24) as f32;
    [
        0.5,
        1.0,
        1.5,
        2.0,
        3.0,
        // Largest non-integer
        P23 - 0.5,
        // Largest odd integers
        P23 - 1.0,
        P23 + 1.0,
        P24 - 1.0,
        // Even integers, all the larger ones are even
        P23,
        P24,
        f32::MAX,
    ]
    .into_iter()
}

/// Returns positive bases at every binade, including subnormals and bases
/// close to one, where `y` needs to be large.
fn test_bases() -> impl Iterator<Item = f32> {
    let near_one = values::binades(-24..=-1, 1000).flat_map(|v| [purify(1.0 + v), purify(1.0 - v)]);
    values::binades(-126..=127, 2000)
        .chain(values::subnormals(10_000))
        .chain(near_one)
}

/// Returns exponents `y` such that `x^y` is spread across the whole range.
///
/// `log2(x^y) = y * log2(x)` takes evenly spread values in `[-160, 140]`
/// (beyond the overflow and underflow limits), tiny values (results close to
/// one), and the values at the limits of the range: the overflow threshold,
/// the smallest normal and subnormal values, and half the smallest subnormal
/// (which is rounded to zero). Each non-integer exponent is also rounded to an
/// integer, which tests negative bases too.
fn result_targets(x: f32) -> Vec<f32> {
    let log2x = fpmath::log2(f64::from(x));
    if log2x == 0.0 {
        return Vec::new();
    }

    let seed = u64::from(x.to_bits());
    let unit = |s: u64| fpmath::scalbn(s as f64, -32);
    let spread = crate::utils::spread(32, 8, seed).map(|s| -160.0 + 300.0 * unit(s));
    let tiny = crate::utils::spread(32, 2, !seed)
        .zip([1.0, -1.0])
        .map(|(s, sign)| sign * (-1.0 - 30.0 * unit(s)).exp2());
    let limits = [128.0, -126.0, -149.0, -150.0];

    let mut ys = Vec::new();
    for t in spread.chain(tiny).chain(limits) {
        let y = (t / log2x) as f32;
        ys.push(y);
        if y.round() != y {
            ys.push(y.round());
        }
    }
    ys
}

/// Returns pairs `(x, y)` with `x = m^d * 2^(d * k)` and `y = n / d`, with odd
/// `m`, so `x^y = m^n * 2^(n * k)` is exact when `m^n` fits in 24 bits, and a
/// tie between two consecutive values when it has 25 bits.
///
/// `k` is chosen so that the results are spread over the whole range, including
/// subnormals (where they are not exact anymore).
fn exact_results() -> Vec<(f32, f32)> {
    let mut cases = Vec::new();
    for d in [1, 2, 4] {
        for n in 1..=25 {
            if n == d || (d != 1 && n % 2 == 0) {
                // `y` is one or an integer
                continue;
            }

            // `m^d` must fit in 24 bits, so `x` is exact, and `m^n` in 25 bits
            let max_m = max_root(d, 24).min(max_root(n, 25));
            if max_m < 3 {
                continue;
            }
            let ms = [3, 5, 7, max_m - 2, max_m]
                .into_iter()
                .chain(
                    crate::utils::spread(32, 32, (d * 100 + n).into()).map(|s| 3 + s % (max_m - 2)),
                )
                .map(|m| m | 1)
                .filter(|&m| m <= max_m);

            for m in ms {
                let md = m.pow(d);
                let mn = m.pow(n);
                let x_bits = 64 - md.leading_zeros() as i32;
                let r_bits = 64 - mn.leading_zeros() as i32;
                let (d, n) = (d as i32, n as i32);

                // `x` is normal and the result is not smaller than half the
                // smallest subnormal or larger than the largest finite value
                let k_lo = ((-126 - x_bits + 1) as f64 / d as f64)
                    .max((-150 - r_bits + 1) as f64 / n as f64)
                    .ceil() as i32;
                let k_hi = ((127 - x_bits + 1) as f64 / d as f64)
                    .min((127 - r_bits + 1) as f64 / n as f64)
                    .floor() as i32;
                if k_lo > k_hi {
                    continue;
                }

                let ks = [k_lo, k_lo + 1, 0, k_hi - 1, k_hi].into_iter().chain(
                    crate::utils::spread(32, 16, m)
                        .map(move |s| k_lo + (s % (k_hi - k_lo + 1) as u64) as i32),
                );
                for k in ks.filter(|k| (k_lo..=k_hi).contains(k)) {
                    let x = md as f32 * exp2i(d * k);
                    let y = n as f32 / d as f32;
                    cases.push((x, y));
                }
            }
        }
    }

    // Powers of two with results at the limits of the range
    for t in [127, 128, -126, -127, -149, -150, -151] {
        for k in (-149..=127).filter(|&k| k != 0 && t % k == 0) {
            cases.push((exp2i(k), (t / k) as f32));
        }
    }

    cases
}

/// Returns the largest `m` such that `m^n < 2^bits`.
fn max_root(n: u32, bits: u32) -> u64 {
    let mut m = 2f64.powf(f64::from(bits) / f64::from(n)) as u64 + 1;
    while m.checked_pow(n).is_none_or(|p| p >= 1 << bits) {
        m -= 1;
    }
    m
}

/// Returns `2^e`, for `e` from -149 to 127.
fn exp2i(e: i32) -> f32 {
    assert!(matches!(e, -149..=127));
    if e >= -126 {
        f32::from_bits(((e + 127) as u32) << 23)
    } else {
        f32::from_bits(1 << (e + 149))
    }
}
