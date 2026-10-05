use super::{MIN_MAX_ERROR, RUG_PREC, check_result, purify, values};

#[test]
fn test_pow() {
    let mut max_error: f64 = 0.0;
    test_pow_with(|x, y| {
        let bigx = rug::Float::with_val(53, x);
        let bigy = rug::Float::with_val(53, y);
        let expected = rug::Float::with_val(RUG_PREC, rug::ops::Pow::pow(&bigx, &bigy));
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
fn check_neg_base(x: f64, y: f64) {
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
    let mut max_error: f64 = 0.0;
    test_powi_with(|x, y| {
        let bigx = rug::Float::with_val(53, x);
        let expected = rug::Float::with_val(RUG_PREC, rug::ops::Pow::pow(&bigx, y));
        let actual = fpmath::powi(x, y);
        assert_result_eq!(fpmath::pow(x, f64::from(y)), actual);
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

fn test_pow_with(mut f: impl FnMut(f64, f64)) {
    // Negative bases are tested by `check_neg_base`.

    // Test special values
    for x in values::specials().chain([1.0]) {
        for y in values::specials().chain(special_exponents()) {
            f(x, y);
            f(x, -y);
        }
    }

    // Test halves and small integers as exponents with bases at many binades
    for x in values::binades((-1022..=1023).step_by(3), 0).chain(values::subnormals(10)) {
        for y in [0.5, 1.0, 1.5, 2.0, 3.0] {
            f(x, y);
            f(x, -y);
        }
    }

    // Test the largest exponents with bases close to one
    for x in values::binades([-1, 0], 0) {
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

fn test_powi_with(mut f: impl FnMut(f64, i32)) {
    // Negative bases are tested by the caller.

    // Test special values
    for x in values::specials().chain([1.0]) {
        for y in [0, 1, 2, 3, i32::MAX - 1, i32::MAX] {
            f(x, y);
            f(x, -y);
        }
        f(x, i32::MIN);
    }

    // Test small exponents with bases at many binades
    for x in values::binades((-1022..=1023).step_by(3), 0).chain(values::subnormals(0)) {
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
fn special_exponents() -> impl Iterator<Item = f64> {
    const P52: f64 = (1u64 << 52) as f64;
    const P53: f64 = (1u64 << 53) as f64;
    [
        0.5,
        1.0,
        1.5,
        2.0,
        3.0,
        // Largest non-integer
        P52 - 0.5,
        // Largest odd integers
        P52 - 1.0,
        P52 + 1.0,
        P53 - 1.0,
        // Even integers, all the larger ones are even
        P52,
        P53,
        f64::MAX,
    ]
    .into_iter()
}

/// Returns positive bases at every binade, including subnormals and bases
/// close to one, where `y` needs to be large.
fn test_bases() -> impl Iterator<Item = f64> {
    let near_one = values::binades(-53..=-1, 20).flat_map(|v| [purify(1.0 + v), purify(1.0 - v)]);
    values::binades(-1022..=1023, 0)
        .chain(values::subnormals(100))
        .chain(near_one)
}

/// Returns exponents `y` such that `x^y` is spread across the whole range.
///
/// `log2(x^y) = y * log2(x)` takes evenly spread values in `[-1100, 1050]`
/// (beyond the overflow and underflow limits), tiny values (results close to
/// one), and the values at the limits of the range: the overflow threshold,
/// the smallest normal and subnormal values, and half the smallest subnormal
/// (which is rounded to zero). Each non-integer exponent is also rounded to an
/// integer, which tests negative bases too.
fn result_targets(x: f64) -> Vec<f64> {
    let log2x = rug::Float::with_val(RUG_PREC, x).log2();
    if log2x.is_zero() {
        return Vec::new();
    }

    let seed = x.to_bits();
    let unit = |s: u64| fpmath::scalbn(s as f64, -32);
    let spread = crate::utils::spread(32, 3, seed).map(|s| -1100.0 + 2150.0 * unit(s));
    let tiny = crate::utils::spread(32, 2, !seed)
        .zip([1.0, -1.0])
        .map(|(s, sign)| sign * (-1.0 - 60.0 * unit(s)).exp2());
    let limits = [1024.0, -1022.0, -1074.0, -1075.0];

    let mut ys = Vec::new();
    for t in spread.chain(tiny).chain(limits) {
        let y = (rug::Float::with_val(RUG_PREC, t) / &log2x).to_f64();
        ys.push(y);
        if y.round() != y {
            ys.push(y.round());
        }
    }
    ys
}

/// Returns pairs `(x, y)` with `x = m^d * 2^(d * k)` and `y = n / d`, with odd
/// `m`, so `x^y = m^n * 2^(n * k)` is exact when `m^n` fits in 53 bits, and a
/// tie between two consecutive values when it has 54 bits.
///
/// `k` is chosen so that the results are spread over the whole range, including
/// subnormals (where they are not exact anymore).
fn exact_results() -> Vec<(f64, f64)> {
    let mut cases = Vec::new();
    for d in [1, 2, 4] {
        for n in 1..=54 {
            if n == d || (d != 1 && n % 2 == 0) {
                // `y` is one or an integer
                continue;
            }

            // `m^d` must fit in 53 bits, so `x` is exact, and `m^n` in 54 bits
            let max_m = max_root(d, 53).min(max_root(n, 54));
            if max_m < 3 {
                continue;
            }
            let ms = [3, 5, 7, max_m - 2, max_m]
                .into_iter()
                .chain(
                    crate::utils::spread(32, 8, (d * 100 + n).into()).map(|s| 3 + s % (max_m - 2)),
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
                let k_lo = ((-1022 - x_bits + 1) as f64 / d as f64)
                    .max((-1075 - r_bits + 1) as f64 / n as f64)
                    .ceil() as i32;
                let k_hi = ((1023 - x_bits + 1) as f64 / d as f64)
                    .min((1023 - r_bits + 1) as f64 / n as f64)
                    .floor() as i32;
                if k_lo > k_hi {
                    continue;
                }

                let ks = [k_lo, k_lo + 1, 0, k_hi - 1, k_hi].into_iter().chain(
                    crate::utils::spread(32, 4, m)
                        .map(move |s| k_lo + (s % (k_hi - k_lo + 1) as u64) as i32),
                );
                for k in ks.filter(|k| (k_lo..=k_hi).contains(k)) {
                    let x = md as f64 * exp2i(d * k);
                    let y = f64::from(n) / f64::from(d);
                    cases.push((x, y));
                }
            }
        }
    }

    // Powers of two with results at the limits of the range
    for t in [1023, 1024, -1022, -1023, -1074, -1075, -1076] {
        for k in (-1074..=1023).filter(|&k| k != 0 && t % k == 0) {
            cases.push((exp2i(k), f64::from(t / k)));
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

/// Returns `2^e`, for `e` from -1074 to 1023.
fn exp2i(e: i32) -> f64 {
    assert!(matches!(e, -1074..=1023));
    if e >= -1022 {
        f64::from_bits(((e + 1023) as u64) << 52)
    } else {
        f64::from_bits(1 << (e + 1074))
    }
}
