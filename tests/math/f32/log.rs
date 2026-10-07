use super::{MIN_MAX_ERROR, check_result, purify, values};

#[test]
fn test_ln() {
    let mut max_error: f32 = 0.0;
    test_log_with(Func::Ln, |x| {
        let expected = fpmath::ln(f64::from(x));
        let actual = fpmath::ln(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_ln_1p() {
    let mut max_error: f32 = 0.0;
    test_ln_1p_with(|x| {
        let expected = fpmath::ln_1p(f64::from(x));
        let actual = fpmath::ln_1p(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_log2() {
    let mut max_error: f32 = 0.0;
    test_log_with(Func::Log2, |x| {
        let expected = fpmath::log2(f64::from(x));
        let actual = fpmath::log2(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_log10() {
    let mut max_error: f32 = 0.0;
    test_log_with(Func::Log10, |x| {
        let expected = fpmath::log10(f64::from(x));
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

impl Func {
    /// Calculates the function in `f64`.
    fn eval(self, x: f64) -> f64 {
        match self {
            Self::Ln => fpmath::ln(x),
            Self::Ln1p => fpmath::ln_1p(x),
            Self::Log2 => fpmath::log2(x),
            Self::Log10 => fpmath::log10(x),
        }
    }

    /// Calculates the inverse of the function in `f64`.
    fn inverse(self, y: f64) -> f64 {
        match self {
            Self::Ln => fpmath::exp(y),
            Self::Ln1p => fpmath::exp_m1(y),
            Self::Log2 => fpmath::exp2(y),
            Self::Log10 => fpmath::exp10(y),
        }
    }
}

fn test_log_with(func: Func, mut f: impl FnMut(f32)) {
    // Test special values
    for x in values::specials() {
        f(x);
        f(-x);
    }

    // Exhaustive test of all subnormal numbers
    for x in values::subnormals_full() {
        f(x);
    }

    // Exhaustive test of all mantissas at the binades around one, which
    // include all the arguments with results close to zero
    for x in values::binades_full(-1..=0) {
        f(x);
    }

    // Test across a wide range of normal numbers
    for x in values::binades(-126..=127, 100_000) {
        f(x);
    }

    // Test the mantissa patterns of negative values (whose results are NaN)
    for x in values::subnormals(0).chain(values::binades(-126..=127, 0)) {
        f(-x);
    }

    // Test the arguments around the centers and boundaries of the
    // subintervals of the table used by the implementation, at each binade
    for x in table_points(-126..=127) {
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
        f(i as f32);
    }
}

fn test_ln_1p_with(mut f: impl FnMut(f32)) {
    // Test special values
    for x in values::specials() {
        f(x);
        f(-x);
    }

    // Test a set of subnormal values, whose results are trivial (`x`)
    for x in values::subnormals(10_000) {
        f(x);
        f(-x);
    }

    // Exhaustive test of all mantissas at the binades around one
    for x in values::binades_full(-1..=0) {
        f(x);
    }

    // Exhaustive test of all the values in (-1, -1/2], whose results are
    // large and negative
    for x in values::binades_full(-1..=-1) {
        f(-x);
    }

    // Test the mantissa patterns of the values from -1 down: -1, whose result
    // is -inf, and the values below it, whose results are NaN (including the
    // closest ones)
    for x in values::binades(0..=127, 0) {
        f(-x);
    }

    // Test across a wide range of normal numbers, positive and in (-1, 0)
    for x in values::binades(-126..=127, 100_000) {
        f(x);
    }
    for x in values::binades(-126..=-1, 100_000) {
        f(-x);
    }

    // Test the arguments around the centers and boundaries of the
    // subintervals of the table used by the implementation (for `1 + x`), at
    // each binade
    for y in table_points(-1..=127) {
        f(purify(y - 1.0));
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
        f(i as f32);
    }
}

/// Number of subintervals in which the table used by the implementation
/// splits each binade.
const TABLE_N: i32 = 128;

/// Returns the arguments around the centers (where the reduced argument is
/// zero) and the boundaries (where it is the largest) of the subintervals of
/// the table used by the implementation, at the binades with exponents in `e`.
fn table_points(e: impl IntoIterator<Item = i32>) -> impl Iterator<Item = f32> {
    e.into_iter().flat_map(|e| {
        (0..2 * TABLE_N).flat_map(move |i| {
            let x = 1.0 + i as f32 / (2 * TABLE_N) as f32;
            values::around(fpmath::scalbn(x, e), 1)
        })
    })
}

/// Returns the arguments around the integer powers of the base (`b^n`, or
/// `e^n - 1` for `ln_1p`) in the range of finite values.
fn base_powers(func: Func) -> impl Iterator<Item = f32> {
    let lowest = match func {
        // The smallest value above -1 is `-1 + 2^-24`
        Func::Ln1p => fpmath::scalbn(1.0, -24) - 1.0,
        _ => fpmath::scalbn(1.0, -149),
    };
    let n_lo = func.eval(lowest).ceil() as i32;
    let n_hi = func.eval(f64::from(f32::MAX)).floor() as i32;
    (n_lo..=n_hi).flat_map(move |n| values::around(func.inverse(f64::from(n)) as f32, 1))
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
/// `values::offset_midpoint`).
fn near_midpoints(func: Func) -> Vec<f32> {
    let eval = |x: f64| func.eval(x);
    let inverse = |y: f64| func.inverse(y);

    let mut xs = Vec::new();

    // Arguments at the binades far from one (with exponents of at least 16 in
    // magnitude), whose results are large
    let large = (-126..=127)
        .filter(|e: &i16| e.abs() >= 16 && (e.is_positive() || !matches!(func, Func::Ln1p)));
    for x0 in values::binades(large, 1000) {
        xs.extend(values::around(
            values::midpoint_inverse(x0, eval, inverse),
            1,
        ));
    }

    if func == Func::Ln1p {
        for x0 in values::binades(-24..=-2, 5000).flat_map(|x| [x, -x]) {
            // g(x) = ln_1p(x) - x, g'(x) = 1 / (1 + x) - 1 = -x / (1 + x)
            if let Some(x) = values::offset_midpoint(x0, eval, |x| -x / (1.0 + x)) {
                xs.extend(values::around(x, 1));
            }
        }
    }

    // Close arguments can give the same midpoint
    xs.sort_by(f32::total_cmp);
    xs.dedup();
    xs
}
