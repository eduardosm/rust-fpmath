use super::{MIN_MAX_ERROR, RUG_PREC, check_result, exponent, purify, values};

#[test]
fn test_ln() {
    let mut max_error: f64 = 0.0;
    test_log_with(Func::Ln, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).ln();
        let actual = fpmath::ln(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_ln_1p() {
    let mut max_error: f64 = 0.0;
    test_ln_1p_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).ln_1p();
        let actual = fpmath::ln_1p(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_log2() {
    let mut max_error: f64 = 0.0;
    test_log_with(Func::Log2, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).log2();
        let actual = fpmath::log2(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_log10() {
    let mut max_error: f64 = 0.0;
    test_log_with(Func::Log10, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).log10();
        let actual = fpmath::log10(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

/// A function of the `log` family.
#[derive(Copy, Clone, PartialEq, Eq)]
enum Func {
    Ln,
    Ln1p,
    Log2,
    Log10,
}

/// Precision used to calculate the test values.
const EXT_PREC: u32 = 256;

impl Func {
    /// Calculates the function with high precision.
    fn eval(self, x: &rug::Float) -> rug::Float {
        let x = rug::Float::with_val(EXT_PREC, x);
        match self {
            Self::Ln => x.ln(),
            Self::Ln1p => x.ln_1p(),
            Self::Log2 => x.log2(),
            Self::Log10 => x.log10(),
        }
    }

    /// Calculates the inverse of the function with high precision.
    fn inverse(self, y: &rug::Float) -> rug::Float {
        let y = rug::Float::with_val(EXT_PREC, y);
        match self {
            Self::Ln => y.exp(),
            Self::Ln1p => y.exp_m1(),
            Self::Log2 => y.exp2(),
            Self::Log10 => y.exp10(),
        }
    }
}

fn test_log_with(func: Func, mut f: impl FnMut(f64)) {
    // Test special values
    for x in values::specials() {
        f(x);
        f(-x);
    }

    // Test a set of subnormal values
    for x in values::subnormals(100_000) {
        f(x);
    }

    // Test a set of normal values at each binade
    for x in values::binades(-1022..=1023, 500) {
        f(x);
    }

    // Test the mantissa patterns of negative values (whose results are NaN)
    for x in values::subnormals(0).chain(values::binades(-1022..=1023, 0)) {
        f(-x);
    }

    // Test more values at the binades around one
    for x in values::binades(-1..=0, 100_000) {
        f(x);
    }

    // Test all the values within 2^16 ULPs of one, whose results are close to
    // zero
    for x in values::around(1.0, 1 << 16) {
        f(x);
    }

    // Test values close to one at larger distances (`1 + v` and `1 - v` with
    // `v` at each binade below one)
    for v in values::binades(-36..=-1, 10_000) {
        f(purify(1.0 + v));
        f(purify(1.0 - v));
    }

    // Test the arguments around the centers and boundaries of the
    // subintervals of the table used by the implementation, at many binades
    for x in table_points((-1022..=1023).filter(|e: &i32| e % 4 == 0 || e.abs() <= 1)) {
        f(x);
    }

    // Test the arguments around integer powers of the base, whose results are
    // close to integers (exact with powers of two for `log2` and with small
    // powers of ten for `log10`)
    for x in base_powers(func) {
        f(x);
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round
    for x in near_midpoints(func) {
        f(x);
    }

    // Test integers
    for i in 1..=10000 {
        f(i as f64);
    }
}

fn test_ln_1p_with(mut f: impl FnMut(f64)) {
    // Test special values
    for x in values::specials() {
        f(x);
        f(-x);
    }

    // Test a set of subnormal values, whose results are trivial (`x`)
    for x in values::subnormals(1000) {
        f(x);
        f(-x);
    }

    // Test a set of normal values at each binade, positive and in (-1, 0)
    for x in values::binades(-1022..=1023, 500) {
        f(x);
    }
    for x in values::binades(-1022..=-1, 500) {
        f(-x);
    }

    // Test more values with small magnitudes, whose results are close to `x`
    for x in values::binades(-60..=-1, 10_000) {
        f(x);
        f(-x);
    }

    // Test more values at the binades around one
    for x in values::binades(-1..=0, 100_000) {
        f(x);
    }

    // Test all the values within 2^16 ULPs above -1, whose results are large
    // and negative
    for x in values::around(-1.0, 1 << 16).filter(|&x| x >= -1.0) {
        f(x);
    }

    // Test the mantissa patterns of the values from -1 down: -1, whose result
    // is -inf, and the values below it, whose results are NaN (including the
    // closest ones)
    for x in values::binades(0..=1023, 0) {
        f(-x);
    }

    // Test values close to -1 at larger distances (`v - 1` with `v` at each
    // binade below 1/2)
    for v in values::binades(-37..=-2, 10_000) {
        f(purify(v - 1.0));
    }

    // Test the arguments around the centers and boundaries of the
    // subintervals of the table used by the implementation (for `1 + x`), at
    // each binade
    for y in table_points((-1..=1023).filter(|e| e % 4 == 0 || *e <= 1)) {
        f(purify(y - 1.0));
    }

    // Test arguments where `1 + x` is not exact and its reduced argument is
    // very close to zero
    for x in table_reciprocals() {
        f(x);
    }

    // Test the arguments around `e^n - 1` for integers `n`, whose results are
    // close to integers
    for x in base_powers(Func::Ln1p) {
        f(x);
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round
    for x in near_midpoints(Func::Ln1p) {
        f(x);
    }

    // Test integers
    for i in 1..=10000 {
        f(i as f64);
    }
}

/// Number of subintervals in which the table used by the implementation
/// splits each binade.
const TABLE_N: i32 = 128;

/// Returns the arguments around the centers (where the reduced argument is
/// zero) and the boundaries (where it is the largest) of the subintervals of
/// the table used by the implementation, at the binades with exponents in `e`.
fn table_points(e: impl IntoIterator<Item = i32>) -> impl Iterator<Item = f64> {
    e.into_iter().flat_map(|e| {
        (0..2 * TABLE_N).flat_map(move |i| {
            let x = 1.0 + f64::from(i) / f64::from(2 * TABLE_N);
            values::around(fpmath::scalbn(x, e), 1)
        })
    })
}

/// Returns the arguments `x` around `2^k / s - 1`, where `s` is a scale factor
/// of the table used by the implementation (`1 / (1 + i / N)` rounded to 24
/// bits), at the binades where `1 + x` can be not exact.
///
/// The reduced argument of `1 + x` is very close to zero, and it can be
/// smaller than the error of `1 + x` (scaled by `2^-k * s`).
fn table_reciprocals() -> impl Iterator<Item = f64> {
    (1..TABLE_N).flat_map(|i| {
        let s = rug::Float::with_val(24, TABLE_N) / (TABLE_N + i);
        let r = rug::Float::with_val(EXT_PREC, 1) / s;
        (-1..=7).chain(53..=63).flat_map(move |k: i32| {
            let y = rug::Float::with_val(EXT_PREC, &r << k).to_f64();
            values::around(purify(y - 1.0), 16)
        })
    })
}

/// Returns the arguments around the integer powers of the base (`b^n`, or
/// `e^n - 1` for `ln_1p`) in the range of finite values.
fn base_powers(func: Func) -> impl Iterator<Item = f64> {
    let lowest = match func {
        // The smallest value above -1 is `-1 + 2^-53`
        Func::Ln1p => (rug::Float::with_val(EXT_PREC, 1) << -53) - 1,
        _ => rug::Float::with_val(EXT_PREC, 1) << -1074,
    };
    let n_lo = func.eval(&lowest).to_f64().ceil() as i32;
    let n_hi = func
        .eval(&rug::Float::with_val(EXT_PREC, f64::MAX))
        .to_f64()
        .floor() as i32;
    (n_lo..=n_hi).flat_map(move |n| {
        values::around(func.inverse(&rug::Float::with_val(EXT_PREC, n)).to_f64(), 1)
    })
}

/// Returns arguments whose results are very close to a midpoint between two
/// consecutive values.
///
/// With large results, the functions are flat (their results change much
/// less than an ULP between consecutive arguments), so the arguments closest
/// to the inverse of a midpoint have such results (see
/// `values::midpoint_inverse`), and they are closer for larger results.
///
/// With small arguments, `ln_1p(x) = x + g(x)` is not flat, but `g` is (see
/// `values::offset_midpoint`). For `ln` close to one, see [`ln_near_one`].
fn near_midpoints(func: Func) -> Vec<f64> {
    let eval = |x: &rug::Float| func.eval(x);
    let inverse = |y: &rug::Float| func.inverse(y);

    let mut xs = Vec::new();

    // Arguments at every other binade far from one (with exponents of at least
    // 64 in magnitude), whose results are large
    let large = (-1022..=1023).filter(|e: &i16| {
        e.abs() >= 64 && e % 2 == 0 && (e.is_positive() || !matches!(func, Func::Ln1p))
    });
    for x0 in values::binades(large, 20) {
        xs.extend(values::around(
            values::midpoint_inverse(x0, eval, inverse),
            1,
        ));
    }

    match func {
        Func::Ln => xs.extend(ln_near_one()),
        Func::Ln1p => {
            for x0 in values::binades(-53..=-2, 500).flat_map(|x| [x, -x]) {
                // g(x) = ln_1p(x) - x, g'(x) = 1 / (1 + x) - 1 = -x / (1 + x)
                let dg = |x: &rug::Float| -> rug::Float {
                    let r: rug::Float = rug::Float::with_val(EXT_PREC, x) / (x.clone() + 1);
                    -r
                };
                if let Some(x) = values::offset_midpoint(x0, eval, dg) {
                    xs.extend(values::around(x, 1));
                }
            }
        }
        Func::Log2 | Func::Log10 => {}
    }

    // Close arguments can give the same midpoint
    xs.sort_by(f64::total_cmp);
    xs.dedup();
    xs
}

/// Returns arguments close to one whose results with `ln` are very close to a
/// midpoint between two consecutive values.
///
/// With `x = 1 + d`, `ln(x) = d - d^2/2 + d^3/3 - ...`. When `d = j * 2^-52`
/// with `j = 2^b * o` and an odd `o` in `[2^b, 2^(b+1))`, `d` is in the binade
/// with exponent `2b - 52`, and `d^2/2` is an odd multiple of half an ULP of
/// `d`. So, when the result is in that binade too, `d - d^2/2` is a midpoint,
/// which the result misses by about `d^3/3`, less than 2^-10 ULPs with
/// `b <= 10`. Below one, the same happens with `x = 1 - d`, `d = j * 2^-53`
/// and an odd `o` in `[2^(b-1), 2^b)`.
fn ln_near_one() -> Vec<f64> {
    let mut xs = Vec::new();
    for b in 0..=10 {
        for o in ((1 << b)..(2 << b)).filter(|o| o % 2 == 1) {
            xs.push(1.0 + fpmath::scalbn(f64::from(o << b), -52));
        }
        if b >= 1 {
            for o in ((1 << (b - 1))..(1 << b)).filter(|o| o % 2 == 1) {
                xs.push(1.0 - fpmath::scalbn(f64::from(o << b), -53));
            }
        }
    }

    // Keep the ones whose results are in the binade of `d`
    xs.retain(|&x| {
        let y = Func::Ln.eval(&rug::Float::with_val(EXT_PREC, x)).to_f64();
        exponent(y) == exponent(x - 1.0)
    });
    xs
}
