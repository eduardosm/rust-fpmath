const SPREAD_SQRT2: u64 = 0x6A09_E667_F3BC_C909; // 2^64 * frac(sqrt(2))

/// Evenly spread values in `[0, 2^bits)` (Weyl sequence based on the golden
/// ratio).
pub(crate) fn spread(bits: u32, n: u32, seed: u64) -> impl Clone + Iterator<Item = u64> {
    const K: u64 = 0x9E37_79B9_7F4A_7C15; // 2^64 / golden ratio
    let start = seed.wrapping_mul(SPREAD_SQRT2);
    (1..=u64::from(n)).map(move |i| start.wrapping_add(i.wrapping_mul(K)) >> (64 - bits))
}

/// Evenly spread pairs of values in `[0, 2^bits)^2` (R2 sequence, the 2D
/// generalization of the golden ratio sequence, based on the plastic number).
///
/// The pairs are spread evenly over the square, unlike the product of two 1D
/// sequences, which repeats each coordinate many times.
pub(crate) fn spread2(bits: u32, n: u32, seed: u64) -> impl Clone + Iterator<Item = (u64, u64)> {
    const K1: u64 = 0xC13F_A9A9_02A6_328F; // 2^64 / plastic number
    const K2: u64 = 0x91E1_0DA5_C79E_7B1C; // 2^64 / plastic number^2
    let start = seed.wrapping_mul(SPREAD_SQRT2);
    (1..=u64::from(n)).map(move |i| {
        (
            start.wrapping_add(i.wrapping_mul(K1)) >> (64 - bits),
            start.wrapping_add(i.wrapping_mul(K2)) >> (64 - bits),
        )
    })
}

/// Returns the integers `m` in `[2^(bits - 1), 2^bits)` such that `m * 2^s`
/// is very close to a multiple of a constant `c`, which `c(prec)` calculates
/// with precision `prec`.
///
/// `m * 2^s - k * c` is small when `m / k` is a good rational approximation of
/// `b = c * 2^-s`, and the best ones are the convergents and semiconvergents
/// of the continued fraction `[a[0]; a[1], a[2], ...]` of `b`. For each
/// convergent `p[i] / q[i]`, this includes the numerators `p[i-1] + t * p[i]`
/// with `t <= a[i+1]` (semiconvergents) and the two largest `t` in range (the
/// closest ones), and the multiples `t * p[i]` with the two smallest `t` in
/// range.
pub(crate) fn near_multiples(c: impl Fn(u32) -> rug::Float, s: i32, bits: u32) -> Vec<u64> {
    let lo = 1u128 << (bits - 1);
    let hi = 1u128 << bits;

    // `b ~= num / den`, precise enough for the continued fraction to be exact
    // up to denominators larger than 2^bits / b
    let prec = 2 * (bits + s.unsigned_abs()) + 128;
    let mut num = (c(prec) << (i64::from(prec) - i64::from(s)) as i32)
        .to_integer()
        .unwrap();
    let mut den = rug::Integer::from(1) << prec;

    let mut ms = Vec::new();
    // p[i-2] and p[i-1]
    let (mut p0, mut p1) = (0u128, 1u128);
    while den != 0 && p1 < hi {
        // a[i]
        let a = rug::Integer::from(&num / &den);
        let rem = rug::Integer::from(&num - &a * &den);
        (num, den) = (den, rem);
        let a = a.to_u128().unwrap_or(u128::MAX);

        // `p1` is zero when `b < 1` (`p[0] = a[0] = 0`)
        if p1 != 0 {
            // Semiconvergents `p0 + t * p1`, with `1 <= t <= a`
            let t_lo = lo.saturating_sub(p0).div_ceil(p1).max(1);
            let t_hi = ((hi - 1).saturating_sub(p0) / p1).min(a);
            for t in [t_hi, t_hi.saturating_sub(1)] {
                if t >= t_lo && t <= t_hi {
                    ms.push((p0 + t * p1) as u64);
                }
            }

            // Multiples `t * p1` (except for `p1 = p[-1] = 1`)
            if p0 != 0 {
                let t_lo = lo.div_ceil(p1);
                for t in [t_lo, t_lo + 1] {
                    if t * p1 < hi {
                        ms.push((t * p1) as u64);
                    }
                }
            }
        }

        let p2 = a.saturating_mul(p1).saturating_add(p0);
        (p0, p1) = (p1, p2);
    }

    ms.sort_unstable();
    ms.dedup();
    ms
}

#[test]
fn test_spread() {
    for (bits, n) in [(8, 50), (32, 1000), (53, 100_000)] {
        let mut v: Vec<u64> = spread(bits, n, 0x1234).collect();
        assert_eq!(v.len(), n as usize);
        assert!(v.iter().all(|&x| x < (1 << bits)));

        // The values are evenly spread: the largest gap between consecutive
        // values (including the one that wraps around) is a small multiple of
        // the average gap
        v.sort_unstable();
        let max_gap = v
            .windows(2)
            .map(|w| w[1] - w[0])
            .chain([v[0] + (1 << bits) - v[v.len() - 1]])
            .max()
            .unwrap();
        assert!(max_gap <= 2 * ((1 << bits) / u64::from(n)) + 2);
    }
}

#[test]
fn test_spread2() {
    // 2^14 pairs in a grid of 2^6 x 2^6 cells, 4 per cell on average
    let bits = 16;
    let n = 1 << 14;
    let mut cells = vec![0u32; 1 << 12];
    for (a, b) in spread2(bits, n, 0x1234) {
        assert!(a < (1 << bits) && b < (1 << bits));
        cells[((a >> (bits - 6)) << 6 | (b >> (bits - 6))) as usize] += 1;
    }

    // The pairs are evenly spread: every cell has a few of them (with random
    // pairs, many cells would be empty)
    assert!(cells.iter().all(|&c| (1..=8).contains(&c)));
}

#[test]
fn test_near_multiples() {
    let constants: [fn(u32) -> rug::Float; 3] = [
        |prec| rug::Float::with_val(prec, rug::float::Constant::Pi) / 2,
        |prec| rug::Float::with_val(prec, 1).exp(),
        |prec| rug::Float::with_val(prec, 2).sqrt(),
    ];
    let bits = 10;
    let (lo, hi) = (1u64 << (bits - 1), 1u64 << bits);

    for c in constants {
        // From values much smaller than `c` to much larger
        for s in -20..=40 {
            let ms = near_multiples(c, s, bits);
            assert!(ms.iter().all(|m| (lo..hi).contains(m)));
            assert!(ms.windows(2).all(|w| w[0] < w[1]));

            // Exhaustive search of the value closest to a (non-zero) multiple
            let prec = 2 * (bits + s.unsigned_abs()) + 128;
            let c = c(prec);
            let dist = |m: u64| {
                let x = rug::Float::with_val(prec, m) << s;
                let k = rug::Float::with_val(prec, &x / &c)
                    .round()
                    .max(&rug::Float::with_val(prec, 1));
                (x - k * &c).abs()
            };
            let closest = (lo..hi)
                .min_by(|&a, &b| dist(a).total_cmp(&dist(b)))
                .unwrap();
            assert!(
                ms.contains(&closest),
                "s = {s}, closest = {closest}, ms = {ms:?}"
            );
        }
    }
}
