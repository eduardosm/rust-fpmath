const MANT_BITS: u32 = 52;
const MANT_MASK: u64 = (1 << MANT_BITS) - 1;
const MIN_EXP: i16 = -1022;
const MAX_EXP: i16 = 1023;
/// Exponent of the smallest positive subnormal value.
const MIN_SUB_EXP: i16 = MIN_EXP - MANT_BITS as i16;

fn mk_normal(m: u64, e: i16, s: bool) -> f64 {
    assert!(m <= MANT_MASK);
    assert!(matches!(e, MIN_EXP..=MAX_EXP));
    let e = u64::from((e + 1023) as u16) << MANT_BITS;
    let s = u64::from(s) << 63;
    f64::from_bits(m | e | s)
}

fn mk_subnormal(m: u64, s: bool) -> f64 {
    assert!(m <= MANT_MASK);
    let s = u64::from(s) << 63;
    f64::from_bits(m | s)
}

/// Returns an iterator over a set of positive special `f64` values,
/// including NaN, infinity, and extreme values.
pub(super) fn specials() -> impl Iterator<Item = f64> {
    [
        f64::NAN,
        f64::INFINITY,
        f64::MAX,
        f64::MIN_POSITIVE,
        0.0,
        mk_subnormal(1, false),
        mk_subnormal(2, false),
    ]
    .into_iter()
}

/// Returns `(1 + m * 2^-52) * 2^e` for `e` from `MIN_SUB_EXP` to `MAX_EXP`. When
/// `e < MIN_EXP`, the value is subnormal and the low bits of the mantissa that
/// do not fit are discarded.
fn mk_value(m: u64, e: i16) -> f64 {
    assert!(m <= MANT_MASK);
    if e >= MIN_EXP {
        mk_normal(m, e, false)
    } else {
        assert!(e >= MIN_SUB_EXP);
        mk_subnormal(((1 << MANT_BITS) | m) >> (MIN_EXP - e), false)
    }
}

/// Mantissa patterns at the ends and the middle of the range.
fn basic_mantissa_patterns() -> [u64; 11] {
    let max = MANT_MASK;
    let alt = max / 3; // 0b0101...
    [
        0,
        1,
        2,
        max,
        max - 1,
        max - 2,
        1 << (MANT_BITS - 1),
        (1 << (MANT_BITS - 1)) - 1,
        (1 << (MANT_BITS - 1)) + 1,
        alt,
        alt << 1,
    ]
}

fn mantissa_patterns() -> impl Iterator<Item = u64> {
    let max = MANT_MASK;
    let mut patterns: Vec<u64> = basic_mantissa_patterns()
        .into_iter()
        .chain((2..=4).flat_map(move |i| {
            // .01b, .11b, .001b, .111b, .0001b, .1111b
            [1 << (MANT_BITS - i), max - ((1 << (MANT_BITS - i)) - 1)]
        }))
        .chain((0..MANT_BITS).flat_map(|i| {
            // .100000...b, .111111...b, .010000...b, .011111...b, ...
            [1 << i, (1 << (i + 1)) - 1]
        }))
        .collect();
    // Some are repeated (such as `1 << (MANT_BITS - 1)`)
    patterns.sort_unstable();
    patterns.dedup();
    patterns.into_iter()
}

/// Returns an iterator over a set of positive subnormal `f64` values,
/// include `n_spread` evenly spread values.
pub(super) fn subnormals(n_spread: u32) -> impl Iterator<Item = f64> {
    mantissa_patterns()
        .filter(|&m| m != 0)
        .chain(spread(n_spread, 0xE43A))
        .map(|m| mk_subnormal(m, false))
}

/// Returns an iterator over a set of positive normal `f64` values at the given
/// binades. Each binade will include a set of mantissa patterns and `n_spread`
/// evenly spread values.
pub(super) fn binades(
    e: impl IntoIterator<Item = i16>,
    n_spread: u32,
) -> impl Iterator<Item = f64> {
    e.into_iter().flat_map(move |e| {
        mantissa_patterns()
            .chain(spread(n_spread, i64::from(e) as u64))
            .map(move |m| mk_normal(m, e, false))
    })
}

/// Evenly spread values in `[0, MANT_MASK]` (Weyl sequence based on the golden
/// ratio).
pub(crate) fn spread(n: u32, seed: u64) -> impl Iterator<Item = u64> {
    crate::utils::spread(MANT_BITS, n, seed)
}

/// Evenly spread pairs of values in `[0, MANT_MASK]^2` (see [`crate::utils::spread2`]).
fn spread2(n: u32, seed: u64) -> impl Iterator<Item = (u64, u64)> {
    crate::utils::spread2(MANT_BITS, n, seed)
}

/// Returns an iterator over `n` values evenly spread over `[lo, hi]`.
pub(super) fn interval(lo: f64, hi: f64, n: u32) -> impl Iterator<Item = f64> {
    assert!(lo <= hi);
    crate::utils::spread(53, n, lo.to_bits() ^ hi.to_bits())
        .map(move |u| lo + (hi - lo) * fpmath::scalbn(u as f64, -53))
}

/// Returns an iterator over `x` and the `n` values before and after it.
pub(super) fn around(x: f64, n: u32) -> impl Iterator<Item = f64> {
    let max = super::ordinal(f64::INFINITY);
    let ord = super::ordinal(x);
    (ord - i64::from(n)..=ord + i64::from(n))
        .filter(move |o| o.abs() <= max)
        .map(super::from_ordinal)
}

/// Returns the positive values in the binades with exponents in `e` that are
/// very close to multiples of a constant `c`, which `c(prec)` calculates with
/// precision `prec` (see [`crate::utils::near_multiples`]).
pub(super) fn near_multiples(
    c: impl Fn(u32) -> rug::Float,
    e: impl IntoIterator<Item = i16>,
) -> impl Iterator<Item = f64> {
    e.into_iter().flat_map(move |e| {
        assert!(e >= MIN_EXP);
        let s = i32::from(e) - MANT_BITS as i32;
        crate::utils::near_multiples(&c, s, MANT_BITS + 1)
            .into_iter()
            .map(move |m| fpmath::scalbn(m as f64, s))
    })
}

/// Precision used by `rug` to calculate test values.
const EXT_PREC: u32 = 256;

/// Returns the two consecutive values `(below, above)` with
/// `below <= y < above` (`above` is infinity when `y >= MAX`).
fn bracket(y: &rug::Float) -> (f64, f64) {
    let below = y.to_f64_round(rug::float::Round::Down);
    (below, super::from_ordinal(super::ordinal(below) + 1))
}

/// Returns the value closest to the inverse (`inv`) of the midpoint between
/// two consecutive values closest to `f(x)`.
///
/// Where `f` is flat (its result changes much less than an ULP between
/// consecutive arguments), the results of `f` at the returned value and its
/// neighbors are very close to that midpoint.
pub(super) fn midpoint_inverse(
    x: f64,
    f: impl Fn(&rug::Float) -> rug::Float,
    inv: impl Fn(&rug::Float) -> rug::Float,
) -> f64 {
    let (below, above) = bracket(&f(&rug::Float::with_val(EXT_PREC, x)));
    let mid = (rug::Float::with_val(EXT_PREC, below) + above) / 2;
    inv(&mid).to_f64()
}

/// Returns the value closest to the one in the binade of `x0` where
/// `g(x) = f(x) - x` is an odd multiple of half an ULP of `x`, if `f` of it is
/// also in that binade, which makes it very close to a midpoint between two
/// consecutive values.
///
/// This works where `f` is not flat but `g` is, such as `f(x) ~= x` with small
/// `x`. `dg` is the derivative of `g`.
pub(super) fn offset_midpoint(
    x0: f64,
    f: impl Fn(&rug::Float) -> rug::Float,
    dg: impl Fn(&rug::Float) -> rug::Float,
) -> Option<f64> {
    assert!(x0.is_normal());
    let g = |x: &rug::Float| f(x) - x;
    let ulp = rug::Float::with_val(EXT_PREC, 1) << (super::exponent(x0) - 52);
    let target = (rug::Float::with_val(EXT_PREC, g(&rug::Float::with_val(EXT_PREC, x0)) / &ulp)
        .floor()
        + 0.5)
        * &ulp;

    // Newton's method, until it converges (it can start far from the
    // solution at the lowest binades, where `g(x0)` is a few ULPs)
    let mut x = rug::Float::with_val(EXT_PREC, x0);
    for _ in 0..100 {
        let dx: rug::Float = (g(&x) - &target) / dg(&x);
        x -= &dx;
        if dx.is_zero() || dx.abs() < rug::Float::with_val(EXT_PREC, x.abs_ref()) >> 128 {
            break;
        }
    }

    let x = x.to_f64();
    let y = f(&rug::Float::with_val(EXT_PREC, x)).to_f64();
    let e = super::exponent(x0);
    (super::exponent(x) == e && super::exponent(y) == e).then_some(x)
}

/// Returns an iterator over the exponent pairs `(e1, e2)` with `e1 - e2 = d`
/// for each `d` in `diffs`, where both are exponents of positive finite values
/// (from -1074 to 1023, see [`binade_pairs`]).
///
/// For each difference, `e1` takes the two lowest and the two highest valid
/// values, the values that make `e1` or `e2` zero, and every `step`-th value
/// in between.
pub(super) fn exp_pairs(
    diffs: impl IntoIterator<Item = i16>,
    step: i16,
) -> impl Iterator<Item = (i16, i16)> {
    assert!(step > 0);
    diffs.into_iter().flat_map(move |d| {
        let lo = MIN_SUB_EXP.max(MIN_SUB_EXP + d);
        let hi = MAX_EXP.min(MAX_EXP + d);
        (lo..=hi)
            .filter(move |&e| e - lo < 2 || hi - e < 2 || e == 0 || e == d || (e - lo) % step == 0)
            .map(move |e| (e, e - d))
    })
}

/// Returns an iterator over pairs of positive finite `f64` values `(x, y)`
/// where `x` is in the binade with exponent `ex` and `y` in the binade with
/// exponent `ey`, for each `(ex, ey)` in `exps`. Exponents below -1022 give
/// subnormal values, with fewer mantissa bits.
///
/// Each exponent pair includes all the combinations of the basic mantissa
/// patterns, each mantissa pattern paired with itself, and `n_spread` mantissa
/// pairs evenly spread over all the combinations.
pub(super) fn binade_pairs(
    exps: impl IntoIterator<Item = (i16, i16)>,
    n_spread: u32,
) -> impl Iterator<Item = (f64, f64)> {
    let basic = basic_mantissa_patterns();
    let pattern_pairs: Vec<(u64, u64)> = basic
        .iter()
        .flat_map(|&mx| basic.iter().map(move |&my| (mx, my)))
        // The basic ones are already paired with themselves
        .chain(
            mantissa_patterns()
                .filter(|m| !basic.contains(m))
                .map(|m| (m, m)),
        )
        .collect();

    exps.into_iter().flat_map(move |(ex, ey)| {
        let seed = (u64::from(ex as u16) << 16) | u64::from(ey as u16);
        pattern_pairs
            .clone()
            .into_iter()
            .chain(spread2(n_spread, seed))
            .map(move |(mx, my)| (mk_value(mx, ex), mk_value(my, ey)))
    })
}

/// Returns the distance from `y` to the closest midpoint between two
/// consecutive values, in ULPs.
#[cfg(test)]
fn midpoint_distance(y: &rug::Float) -> f64 {
    let (below, above) = bracket(y);
    let frac =
        rug::Float::with_val(EXT_PREC, y - below) / (rug::Float::with_val(EXT_PREC, above) - below);
    (frac.to_f64() - 0.5).abs()
}

#[test]
fn test_mk_value() {
    for e in MIN_SUB_EXP..=MAX_EXP {
        for m in mantissa_patterns() {
            let x = mk_value(m, e);
            // In the binade with exponent `e`
            assert!(x >= fpmath::scalbn(1.0, e.into()));
            assert!(x < fpmath::scalbn(1.0, i32::from(e) + 1));
            // `(1 + m * 2^-52) * 2^e` with the bits that do not fit discarded
            let exact = rug::Float::with_val(64, (1 << MANT_BITS) | m) << (i32::from(e) - 52);
            let ulp = fpmath::scalbn(1.0, (i32::from(e) - 52).max(-1074));
            assert!(x <= exact && exact < rug::Float::with_val(64, x) + ulp);
        }
    }
}

#[test]
fn test_subnormals_not_zero() {
    assert!(subnormals(20_000_000).all(|x| x != 0.0));
}

#[test]
fn test_exp_pairs() {
    for d in [-2097, -2096, -1100, -1, 0, 1, 30, 2096, 2097] {
        // Valid exponents with `e1 - e2 = d`
        let lo = MIN_SUB_EXP.max(MIN_SUB_EXP + d);
        let hi = MAX_EXP.min(MAX_EXP + d);
        for step in [1, 7, 100, 5000] {
            let e1s: Vec<i16> = exp_pairs([d], step)
                .map(|(e1, e2)| {
                    assert_eq!(e1 - e2, d);
                    e1
                })
                .collect();
            assert!(e1s.windows(2).all(|w| w[0] < w[1]));
            assert!(e1s.iter().all(|e| (lo..=hi).contains(e)));

            // The edges, the ones that make `e1` or `e2` zero, and every
            // `step`-th one
            let expected = [lo, lo + 1, hi - 1, hi, 0, d]
                .into_iter()
                .chain((lo..=hi).step_by(step as usize))
                .filter(|e| (lo..=hi).contains(e));
            for e in expected {
                assert!(e1s.contains(&e), "d = {d}, step = {step}, e = {e}");
            }
            if step == 1 {
                assert_eq!(e1s.len(), (hi - lo + 1) as usize);
            }
        }
    }
}

#[test]
fn test_binade_pairs() {
    let n_patterns = mantissa_patterns().count();
    let in_binade = |x: f64, e: i16| {
        x >= fpmath::scalbn(1.0, e.into()) && x < fpmath::scalbn(1.0, i32::from(e) + 1)
    };
    for (ex, ey) in [
        (MIN_SUB_EXP, MIN_SUB_EXP),
        (MIN_SUB_EXP, MAX_EXP),
        (-1023, -1022),
        (0, 0),
        (5, -3),
        (MAX_EXP, MAX_EXP),
    ] {
        let pairs: Vec<(f64, f64)> = binade_pairs([(ex, ey)], 100).collect();
        // All the combinations of the 11 basic patterns, the other patterns
        // paired with themselves, and the spread pairs
        assert_eq!(pairs.len(), 11 * 11 + (n_patterns - 11) + 100);
        assert!(
            pairs
                .iter()
                .all(|&(x, y)| in_binade(x, ex) && in_binade(y, ey))
        );
    }
}

#[test]
fn test_interval() {
    for (lo, hi) in [(0.0, 1.0), (-745.5, 709.5), (1e300, 1.5e300), (2.0, 2.0)] {
        let n = 1000;
        let mut xs: Vec<f64> = interval(lo, hi, n).collect();
        assert_eq!(xs.len(), n as usize);
        assert!(xs.iter().all(|x| (lo..=hi).contains(x)));

        // Evenly spread
        xs.sort_by(f64::total_cmp);
        let max_gap = [xs[0] - lo, hi - xs[xs.len() - 1]]
            .into_iter()
            .chain(xs.windows(2).map(|w| w[1] - w[0]))
            .fold(0.0, f64::max);
        assert!(max_gap <= 2.0 * (hi - lo) / f64::from(n));
    }
}

#[test]
fn test_around() {
    let collect = |x: f64, n: u32| around(x, n).collect::<Vec<f64>>();
    let one_below = 1.0 - f64::EPSILON / 2.0;
    assert_eq!(
        collect(1.0, 2),
        [
            one_below - f64::EPSILON / 2.0,
            one_below,
            1.0,
            1.0 + f64::EPSILON,
            1.0 + 2.0 * f64::EPSILON
        ],
    );
    // Across zero (only one zero)
    let min_sub = fpmath::scalbn(1.0, -1074);
    assert_eq!(collect(0.0, 1), [-min_sub, 0.0, min_sub]);
    assert_eq!(collect(-0.0, 1), [-min_sub, 0.0, min_sub]);
    // Not beyond infinity
    let below_max = super::from_ordinal(super::ordinal(f64::MAX) - 1);
    let below_below_max = super::from_ordinal(super::ordinal(below_max) - 1);
    assert_eq!(
        collect(f64::MAX, 2),
        [below_below_max, below_max, f64::MAX, f64::INFINITY]
    );
    assert_eq!(
        collect(f64::NEG_INFINITY, 1),
        [f64::NEG_INFINITY, -f64::MAX]
    );
}

#[test]
fn test_near_multiples() {
    let half_pi = |prec| rug::Float::with_val(prec, rug::float::Constant::Pi) / 2;
    for e in [0, 1, 52, 100, 849, 1023] {
        let xs: Vec<f64> = near_multiples(half_pi, [e]).collect();
        assert!(!xs.is_empty());
        assert!(xs.iter().all(|&x| {
            x >= fpmath::scalbn(1.0, e.into()) && x < fpmath::scalbn(1.0, i32::from(e) + 1)
        }));
    }

    // The closest value to a multiple of π/2 (at a distance of about
    // 2^-60.89), in the binade with exponent 849
    let worst = fpmath::scalbn(6381956970095103.0, 797);
    assert!(near_multiples(half_pi, [849]).any(|x| x == worst));
}

#[test]
fn test_midpoint_inverse() {
    // The cosine is flat close to zero, its result changes less than
    // `2 * x^2` ULPs between consecutive arguments
    let cos = |x: &rug::Float| rug::Float::with_val(EXT_PREC, x).cos();
    let acos = |y: &rug::Float| rug::Float::with_val(EXT_PREC, y).acos();
    for x0 in binades(-26..=-10, 100) {
        let x = midpoint_inverse(x0, cos, acos);
        let step = 2.0 * x * x;
        let d = midpoint_distance(&cos(&rug::Float::with_val(EXT_PREC, x)));
        assert!(d <= step / 2.0, "x0 = {x0:e}, x = {x:e}, d = {d:e}");
        // The neighbors are on both sides of the midpoint
        for x in [around(x, 1).next().unwrap(), around(x, 1).last().unwrap()] {
            let d = midpoint_distance(&cos(&rug::Float::with_val(EXT_PREC, x)));
            assert!(d <= 1.5 * step, "x0 = {x0:e}, x = {x:e}, d = {d:e}");
        }
    }
}

#[test]
fn test_offset_midpoint() {
    // exp_m1(x) = x + g(x), with g'(x) = exp_m1(x) ~= x
    let exp_m1 = |x: &rug::Float| rug::Float::with_val(EXT_PREC, x).exp_m1();
    for e in -50..=-4 {
        // Away from the edges of the binade, the solution and the result are
        // in the binade
        let lo = fpmath::scalbn(1.25, e);
        let hi = fpmath::scalbn(1.75, e);
        for x0 in interval(lo, hi, 20).flat_map(|x| [x, -x]) {
            let x = offset_midpoint(x0, exp_m1, exp_m1).unwrap();
            assert_eq!(super::exponent(x), super::exponent(x0));
            assert_eq!(x.is_sign_negative(), x0.is_sign_negative());
            // `g` changes by about `x` ULPs between consecutive arguments, so
            // the result is closer than that to a midpoint
            let d = midpoint_distance(&exp_m1(&rug::Float::with_val(EXT_PREC, x)));
            assert!(d < x.abs(), "x0 = {x0:e}, x = {x:e}, d = {d:e}");
        }
    }

    // The result is in the binade below
    assert_eq!(offset_midpoint(-0.0625, exp_m1, exp_m1), None);
    // The result is in the binade above
    let below_eighth = super::from_ordinal(super::ordinal(0.125) - 1);
    assert_eq!(offset_midpoint(below_eighth, exp_m1, exp_m1), None);
}
