use super::values;

#[test]
fn test_round() {
    test_round_with(|arg| {
        let actual = fpmath::round(arg);
        let expected = rug::Float::with_val(128, arg).round().to_f64();
        assert_result_eq!(expected, actual);
    });
}

#[test]
fn test_floor() {
    test_round_with(|arg| {
        let actual = fpmath::floor(arg);
        let expected = rug::Float::with_val(128, arg).floor().to_f64();
        assert_result_eq!(expected, actual);
    });
}

#[test]
fn test_ceil() {
    test_round_with(|arg| {
        let actual = fpmath::ceil(arg);
        let expected = rug::Float::with_val(128, arg).ceil().to_f64();
        assert_result_eq!(expected, actual);
    });
}

#[test]
fn test_trunc() {
    test_round_with(|arg| {
        let actual = fpmath::trunc(arg);
        let expected = rug::Float::with_val(128, arg).trunc().to_f64();
        assert_result_eq!(expected, actual);
    });
}

fn test_round_with(f: fn(f64)) {
    // Test special values
    for x in values::specials() {
        f(x);
        f(-x);
    }

    // Test a set of subnormal values, whose results are zero or one (with the
    // sign of the argument)
    for x in values::subnormals(10_000) {
        f(x);
        f(-x);
    }

    // Test across a set of normal values at each binade
    for x in values::binades(-1022..=1023, 1000) {
        f(x);
        f(-x);
    }

    // Test the values with special integer and fractional parts at each binade
    // where they have both
    for x in split_values() {
        f(x);
        f(-x);
    }
}

/// Returns positive values with special combinations of integer and
/// fractional parts at each binade from one to 2^52 (smaller values are
/// rounded to zero or one, and larger values are integers).
///
/// The integer parts are the smallest and the largest ones (when the
/// result is rounded up, it is in the next binade), and the fractional parts
/// are zero, the smallest and largest ones, and around one half (the ties of
/// `round`). They are combined with each other and with evenly spread integer
/// and fractional parts.
fn split_values() -> Vec<f64> {
    let mut xs = Vec::new();
    for e in 0..=52u32 {
        // Number of fractional bits
        let f = 52 - e;
        let mask = |bits: u32| (1u64 << bits) - 1;
        let spread = |bits: u32, seed: u64| {
            (bits != 0)
                .then(|| crate::utils::spread(bits, 100, seed))
                .into_iter()
                .flatten()
        };

        // Integer parts below the implicit bit
        let ints: Vec<u64> = [0, 1, mask(e) >> 1, mask(e).saturating_sub(1), mask(e)]
            .into_iter()
            .chain(spread(e, 1))
            .collect();
        // Fractional parts
        let half = (1u64 << f) >> 1;
        let fracs: Vec<u64> = [
            0,
            1,
            half.saturating_sub(1),
            half,
            half + 1,
            mask(f).saturating_sub(1),
            mask(f),
        ]
        .into_iter()
        .chain(spread(f, 2))
        .collect();

        for &int in &ints {
            for &frac in &fracs {
                // Values in range (some combinations are repeated with few bits)
                if int <= mask(e) && frac <= mask(f) {
                    let m = (1 << 52) | (int << f) | frac;
                    xs.push(fpmath::scalbn(m as f64, e as i32 - 52));
                }
            }
        }
    }
    xs
}
