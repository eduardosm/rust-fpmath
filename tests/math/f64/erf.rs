use super::{MIN_MAX_ERROR, RUG_PREC, check_result, values};

#[test]
fn test_erf() {
    let mut max_error: f64 = 0.0;
    test_with(Func::Erf, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).erf();
        let actual = fpmath::erf(x);
        assert_result_eq!(fpmath::erf(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_erfc() {
    let mut max_error: f64 = 0.0;
    test_with(Func::Erfc, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).erfc();
        let actual = fpmath::erfc(x);
        check_result(x, actual, expected, &mut max_error);

        let expected = rug::Float::with_val(RUG_PREC, -x).erfc();
        let actual = fpmath::erfc(-x);
        check_result(-x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_erfcx() {
    let mut max_error: f64 = 0.0;
    test_erfcx_with(|x| {
        let expected = rug_erfcx(&rug::Float::with_val(53, x), RUG_PREC);
        let actual = fpmath::erfcx(x);
        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_rug_erfcx() {
    // The continued fraction that `rug_erfcx` uses with large arguments must
    // match the definition of `erfcx`
    for x in [4.0, 4.5, 6.0, 9.0, 15.0, 27.0, 100.0] {
        let x = rug::Float::with_val(53, x);
        let exact = rug::Float::with_val(2 * EXT_PREC, x.square_ref()).exp()
            * rug::Float::with_val(2 * EXT_PREC, &x).erfc();
        let error = (rug_erfcx(&x, EXT_PREC) - &exact) / exact;
        assert!(
            error.abs() < rug::Float::with_val(53, 1) >> EXT_PREC,
            "x = {x}"
        );
    }
}

/// An error function.
#[derive(Copy, Clone, PartialEq, Eq)]
enum Func {
    Erf,
    Erfc,
}

/// Precision used to calculate the test values.
const EXT_PREC: u32 = 256;

impl Func {
    /// Calculates the function with high precision.
    fn eval(self, x: &rug::Float) -> rug::Float {
        let x = rug::Float::with_val(EXT_PREC, x);
        match self {
            Self::Erf => x.erf(),
            Self::Erfc => x.erfc(),
        }
    }

    /// Calculates the derivative of the function with high precision, with
    /// `erf'(x) = 2 / sqrt(π) * exp(-x^2)` and `erfc'(x) = -erf'(x)`.
    fn derivative(self, x: &rug::Float) -> rug::Float {
        let pi = rug::Float::with_val(EXT_PREC, rug::float::Constant::Pi);
        let d = 2 * (-rug::Float::with_val(EXT_PREC, x.square_ref())).exp() / pi.sqrt();
        match self {
            Self::Erf => d,
            Self::Erfc => -d,
        }
    }

    /// Solves `func(x) = y` with Newton's method, starting from `x0`.
    fn solve(self, y: &rug::Float, x0: f64) -> rug::Float {
        let mut x = rug::Float::with_val(EXT_PREC, x0);
        for _ in 0..100 {
            let dx: rug::Float = (self.eval(&x) - y) / self.derivative(&x);
            x -= &dx;
            if dx.is_zero() || dx.abs() < rug::Float::with_val(EXT_PREC, x.abs_ref()) >> 128 {
                break;
            }
        }
        x
    }
}

/// Arguments from which `erf` and `erfc` (of the negated argument) are
/// calculated from a less precise `erfc(|x|)`, which is very small.
const PRECISE_MAX: f64 = 3.5;

/// Arguments from which `erf` rounds to one and `erfc` (of the negated
/// argument) rounds to two.
const SATURATION_MAX: f64 = 6.0;

/// Arguments from which `erfc` rounds to zero.
const ERFC_ZERO_MAX: f64 = 27.25;

/// Arguments from which the results of `erfc` are below 2^-1044: they are
/// subnormal, with at most 30 significant bits, while the relative change of
/// `erfc` between consecutive arguments is about 2^-42 there, so the results
/// change much less than an ULP.
const ERFC_FLAT_MIN: f64 = 26.83;

/// `rug` is very slow calculating `erfc` with arguments in this range (it
/// calculates `1 - erf(x)` with enough precision for the cancellation), so
/// fewer of them are tested.
const ERFC_SLOW_RANGE: std::ops::Range<f64> = 8.0..16.0;

fn test_with(func: Func, mut f: impl FnMut(f64)) {
    // Negative arguments are tested by the callers.

    // The largest argument with a non-trivial result
    let x_max = match func {
        Func::Erf => SATURATION_MAX,
        Func::Erfc => ERFC_ZERO_MAX,
    };

    // Test special values
    for x in values::specials() {
        f(x);
    }

    // Test a set of subnormal values, whose results are subnormal with `erf`
    for x in values::subnormals(10_000) {
        f(x);
    }

    // Test a set of normal values at each binade. The results of the large
    // ones are trivial (zero, one or two), and so are the results of `erfc`
    // with the tiny ones (one). The values closest to each power of two are
    // included, which covers the thresholds where the calculation changes
    // with small arguments (such as 2^-960, 2^-55 and 2^-35).
    for x in values::binades(-1022..=1023, 100) {
        f(x);
    }

    // Test more values with small magnitudes, whose results are calculated
    // with the same polynomial
    for x in values::binades(-60..=-2, 10_000) {
        f(x);
    }

    // Test values evenly spread over the whole range of non-trivial results
    // (fewer of them in the slow range of `erfc`)
    let intervals: &[(f64, f64, u32)] = match func {
        Func::Erf => &[(0.0, x_max + 0.25, 2_000_000)],
        Func::Erfc => &[
            (0.0, ERFC_SLOW_RANGE.start, 400_000),
            (ERFC_SLOW_RANGE.start, ERFC_SLOW_RANGE.end, 40_000),
            (ERFC_SLOW_RANGE.end, x_max + 0.25, 400_000),
        ],
    };
    for &(lo, hi, n) in intervals {
        for x in values::interval(lo, hi, n) {
            f(x);
        }
    }

    // Test the arguments around the limits of the range
    for x in limits() {
        for x in values::around(x, 10_000) {
            f(x);
        }
    }

    // Test the arguments around the limits of the intervals of the
    // approximations of `erfcx`
    for x in interval_limits() {
        for x in values::around(x, 100) {
            f(x);
        }
    }

    // Test the arguments around the boundaries of the argument reduction of
    // the exponential (only some of them in the slow range of `erfc`)
    for (i, x) in reduction_boundaries(x_max).into_iter().enumerate() {
        if func == Func::Erfc && ERFC_SLOW_RANGE.contains(&x) && i % 16 != 0 {
            continue;
        }
        for x in values::around(x, 2) {
            f(x);
        }
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round
    for x in near_midpoints(func) {
        f(x);
    }
}

/// Returns the arguments at the limits of the range: where the functions
/// change the evaluation method (0.5, `PRECISE_MAX` and the saturation limits),
/// where `erfc` rounds to one with tiny arguments (positive and negative),
/// where `erf` rounds to one and `erfc` rounds to two (with the negated
/// argument), and where `erfc` becomes subnormal and rounds to zero.
fn limits() -> [f64; 10] {
    let erfc_inv = |e: i32, x0: f64| Func::Erfc.solve(&exp2i(e), x0).to_f64();
    [
        0.5,
        PRECISE_MAX,
        SATURATION_MAX,
        ERFC_ZERO_MAX,
        erf_inv_tiny(-54),
        erf_inv_tiny(-53),
        erfc_inv(-54, 5.9),
        erfc_inv(-53, 5.8),
        erfc_inv(-1022, 26.5),
        erfc_inv(-1075, 27.2),
    ]
}

/// Returns `2^e` with high precision.
fn exp2i(e: i32) -> rug::Float {
    rug::Float::with_val(EXT_PREC, 1) << e
}

/// Returns the tiny argument where `erf` is `2^e`
/// (`erf(x) ~= x * 2 / sqrt(π)`).
fn erf_inv_tiny(e: i32) -> f64 {
    let pi = rug::Float::with_val(EXT_PREC, rug::float::Constant::Pi);
    (pi.sqrt() * exp2i(e - 1)).to_f64()
}

/// Returns the limits of the intervals of the approximations of `erfcx`:
/// each binade from 0.5 to 32 is split in four intervals of the same length.
fn interval_limits() -> Vec<f64> {
    (-1..=4)
        .flat_map(|e| (0..4).map(move |j| fpmath::scalbn(1.0 + 0.25 * f64::from(j), e)))
        .collect()
}

/// Number of parts in which the argument reduction of `exp` splits each power
/// of two: the arguments are reduced to `m * ln(2) / N + r`, with an integer
/// `m` and `|r| <= ln(2) / (2 * N)`.
const REDUCTION_N: u32 = 128;

/// Returns the arguments up to `x_max` where `m` changes in the argument
/// reduction of `exp(-x^2)`, `sqrt((m + 1/2) * ln(2) / N)`.
fn reduction_boundaries(x_max: f64) -> Vec<f64> {
    let ln_2 = rug::Float::with_val(EXT_PREC, 2).ln();
    let m_max = (x_max * x_max / std::f64::consts::LN_2 * f64::from(REDUCTION_N)) as u32;
    (0..=m_max)
        .map(|m| {
            let x2: rug::Float = rug::Float::with_val(EXT_PREC, 2 * m + 1) * &ln_2;
            (x2 / (2 * REDUCTION_N)).sqrt().to_f64()
        })
        .collect()
}

/// Returns arguments whose results are very close to a midpoint between two
/// consecutive values.
///
/// `erf` is flat with large arguments, and so is `erfc` with large negative
/// arguments, with small ones and from `ERFC_FLAT_MIN` (their results change
/// much less than an ULP between consecutive arguments), so the arguments
/// closest to the inverse of a midpoint have such results (see
/// `values::midpoint_inverse`). The negative arguments of `erfc` are returned
/// negated.
///
/// Otherwise, such results are searched among the arguments close to some
/// values (see `values::midpoint_scan`).
fn near_midpoints(func: Func) -> Vec<f64> {
    let eval = |x: &rug::Float| func.eval(x);
    // `erfc(-x)`
    let eval_neg = |x: &rug::Float| func.eval(&rug::Float::with_val(EXT_PREC, -x));

    let mut xs = Vec::new();
    match func {
        Func::Erf => {
            for x0 in values::binades(-40..=1, 1000) {
                xs.push(values::midpoint_scan(x0, eval, 1 << 12));
            }
            for x0 in values::interval(2.5, 5.93, 20_000) {
                let x = values::midpoint_inverse(x0, eval, |y| func.solve(y, x0));
                xs.extend(values::around(x, 1));
            }
        }
        Func::Erfc => {
            for x0 in values::binades(-54..=-6, 500) {
                let x = values::midpoint_inverse(x0, eval, |y| func.solve(y, x0));
                xs.extend(values::around(x, 1));
                let x = values::midpoint_inverse(x0, eval_neg, |y| -func.solve(y, -x0));
                xs.extend(values::around(x, 1));
            }
            for x0 in values::binades(-5..=1, 1000) {
                xs.push(values::midpoint_scan(x0, eval, 1 << 12));
                xs.push(values::midpoint_scan(x0, eval_neg, 1 << 12));
            }
            for x0 in values::interval(2.5, 5.87, 20_000) {
                let x = values::midpoint_inverse(x0, eval_neg, |y| -func.solve(y, -x0));
                xs.extend(values::around(x, 1));
            }
            for x0 in values::interval(4.0, ERFC_FLAT_MIN, 1500) {
                xs.push(values::midpoint_scan(x0, eval, 1 << 12));
            }
            for x0 in values::interval(ERFC_FLAT_MIN, ERFC_ZERO_MAX, 500) {
                let x = values::midpoint_inverse(x0, eval, |y| func.solve(y, x0));
                xs.extend(values::around(x, 1));
            }
        }
    }

    // Close arguments can give the same midpoint
    xs.sort_by(f64::total_cmp);
    xs.dedup();
    xs
}

/// Calculates `erfcx(x) = exp(x^2) * erfc(x)` with a precision of `prec`
/// bits.
///
/// `rug` does not provide it, and `erfc(x)` underflows when `x` is very
/// large, so its continued fraction is used from 4 (see `rug_erfcx_cf`), and
/// `erfcx(x) = 2 * exp(x^2) - erfcx(-x)` up to -4.
fn rug_erfcx(x: &rug::Float, prec: u32) -> rug::Float {
    // Guard bits for the rounding errors
    let prec = prec + 64;
    // exact
    let x2 = || rug::Float::with_val(2 * x.prec(), x.square_ref());
    if *x >= 4 {
        rug_erfcx_cf(x, prec)
    } else if *x <= -4 {
        let pos = rug_erfcx_cf(&rug::Float::with_val(x.prec(), -x), prec);
        2 * rug::Float::with_val(prec, x2().exp_ref()) - pos
    } else {
        rug::Float::with_val(prec, x2().exp_ref()) * rug::Float::with_val(prec, x).erfc()
    }
}

/// Calculates `erfcx(x)` for `x >= 4` with the continued fraction
///
/// `erfcx(x) = 1 / (sqrt(π) * (x + (1/2) / (x + (2/2) / (x + (3/2) / ...))))`
///
/// whose number of terms is doubled until the result does not change with
/// precision `prec` (it converges faster the larger `x` is).
fn rug_erfcx_cf(x: &rug::Float, prec: u32) -> rug::Float {
    if x.is_infinite() {
        return rug::Float::new(prec);
    }
    // The denominator, with `n` terms
    let eval = |n: u32| {
        let mut t = rug::Float::with_val(prec, x);
        for k in (1..=n).rev() {
            t = rug::Float::with_val(prec, k) / (2 * t) + x;
        }
        t
    };
    let mut n = 16;
    let mut t = eval(n);
    loop {
        n *= 2;
        let next = eval(n);
        let diff = rug::Float::with_val(prec, &t - &next).abs();
        let converged = diff <= rug::Float::with_val(prec, &next >> (prec - 16));
        t = next;
        if converged {
            break;
        }
    }
    1 / (rug::Float::with_val(prec, rug::float::Constant::Pi).sqrt() * t)
}

/// Solves `erfcx(x) = y` with Newton's method, starting from `x0`, with
/// `erfcx'(x) = 2 * x * erfcx(x) - 2 / sqrt(π)`.
fn erfcx_solve(y: &rug::Float, x0: f64) -> rug::Float {
    let pi = rug::Float::with_val(EXT_PREC, rug::float::Constant::Pi);
    let k: rug::Float = 2 / pi.sqrt();
    let mut x = rug::Float::with_val(EXT_PREC, x0);
    for _ in 0..100 {
        let fx = rug_erfcx(&x, EXT_PREC);
        let d: rug::Float = rug::Float::with_val(EXT_PREC, &x * &fx) * 2 - &k;
        let dx: rug::Float = (fx - y) / d;
        x -= &dx;
        if dx.is_zero() || dx.abs() < rug::Float::with_val(EXT_PREC, x.abs_ref()) >> 128 {
            break;
        }
    }
    x
}

/// Arguments from which `erfcx` is calculated with its asymptotic series.
const ERFCX_TBL_MAX: f64 = 28.0;

/// Arguments from which `erfcx` (of the negated argument) is calculated as
/// `2 * exp(x^2)`.
const ERFCX_EXP_MAX: f64 = 7.0;

/// Arguments from which `erfcx` (of the negated argument) overflows.
const ERFCX_INF_MAX: f64 = 26.75;

fn test_erfcx_with(mut f: impl FnMut(f64)) {
    // Both signs of a value
    let both = |x: f64| [x, -x];

    // Test special values
    for x in values::specials().flat_map(both) {
        f(x);
    }

    // Test a set of subnormal values, whose results are one
    for x in values::subnormals(10_000).flat_map(both) {
        f(x);
    }

    // Test a set of normal values at each binade. The results of the tiny
    // ones are one, the results of the large negative ones overflow, and the
    // results of the largest positive ones are subnormal. The values closest
    // to each power of two are included, which covers the thresholds where
    // the calculation changes (such as 2^-55 and 2^35).
    for x in values::binades(-1022..=1023, 100).flat_map(both) {
        f(x);
    }

    // Test more values with small magnitudes
    for x in values::binades(-60..=-2, 5000).flat_map(both) {
        f(x);
    }

    // Test values evenly spread up to the overflow with negative arguments
    // and up to where the asymptotic series is used
    for x in values::interval(0.0, 0.5, 200_000).flat_map(both) {
        f(x);
    }
    for x in values::interval(0.5, ERFCX_TBL_MAX + 0.25, 600_000) {
        f(x);
    }
    for x in values::interval(0.5, ERFCX_INF_MAX, 600_000) {
        f(-x);
    }

    // Test more values where the polynomial of the asymptotic series is not
    // constant
    for x in values::binades(4..=35, 10_000) {
        f(x);
    }

    // Test the arguments around the limits of the range
    for x in erfcx_limits() {
        for x in values::around(x, 10_000) {
            f(x);
        }
    }

    // Test the arguments around the limits of the intervals of the
    // approximations of `erfcx`
    for x in interval_limits().into_iter().flat_map(both) {
        for x in values::around(x, 100) {
            f(x);
        }
    }

    // Test the arguments around the boundaries of the argument reduction of
    // the exponential, which is only used with negative arguments and with
    // small ones (only some of them where the result is just `2 * exp(x^2)`)
    for (i, x) in reduction_boundaries(ERFCX_INF_MAX).into_iter().enumerate() {
        if x >= ERFCX_EXP_MAX && i % 16 != 0 {
            continue;
        }
        for x in values::around(x, 2) {
            f(-x);
            if x < 0.5 {
                f(x);
            }
        }
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round
    for x in erfcx_near_midpoints() {
        f(x);
    }
}

/// Returns the arguments at the limits of the range of `erfcx`: where it
/// changes the evaluation method (also with negative arguments), where it
/// rounds to one with tiny arguments (positive and negative), where it
/// becomes subnormal and where it overflows.
fn erfcx_limits() -> [f64; 13] {
    let pi = rug::Float::with_val(EXT_PREC, rug::float::Constant::Pi);
    // erfcx(x) ~= 1 / (sqrt(π) * x) = 2^-1022
    let subnormal = (exp2i(1022) / pi.sqrt()).to_f64();
    // 2 * exp(x^2) = 2^1024
    let overflow = (rug::Float::with_val(EXT_PREC, 2).ln() * 1023u32)
        .sqrt()
        .to_f64();
    [
        0.5,
        ERFCX_TBL_MAX,
        fpmath::scalbn(1.0, 35),
        subnormal,
        // erfcx(x) ~= 1 - erf(x) with tiny x
        erf_inv_tiny(-54),
        erf_inv_tiny(-53),
        -erf_inv_tiny(-54),
        -erf_inv_tiny(-53),
        -0.5,
        -PRECISE_MAX,
        -ERFCX_EXP_MAX,
        -overflow,
        -ERFCX_INF_MAX,
    ]
}

/// Returns arguments whose results with `erfcx` are very close to a midpoint
/// between two consecutive values.
///
/// `erfcx` is flat with small arguments (its result changes much less than
/// an ULP between consecutive arguments), so the arguments closest to the
/// inverse of a midpoint have such results (see `values::midpoint_inverse`).
///
/// Otherwise, such results are searched among the arguments close to some
/// values (see `values::midpoint_scan`), up to the largest ones, whose results
/// are subnormal.
fn erfcx_near_midpoints() -> Vec<f64> {
    let eval = |x: &rug::Float| rug_erfcx(x, EXT_PREC);

    let mut xs = Vec::new();
    for x0 in values::binades(-54..=-6, 300).flat_map(|x| [x, -x]) {
        let x = values::midpoint_inverse(x0, eval, |y| erfcx_solve(y, x0));
        xs.extend(values::around(x, 1));
    }

    let small = values::binades(-5..=-2, 500).flat_map(|x| [x, -x]);
    let table = values::interval(0.5, ERFCX_TBL_MAX, 1500);
    let negative = values::interval(0.5, 26.6, 1500).map(|x| -x);
    let large = values::binades(4..=35, 30).chain(values::binades((36..=1023).step_by(8), 10));
    let subnormal = values::binades(1021..=1023, 200);
    for x0 in small
        .chain(table)
        .chain(negative)
        .chain(large)
        .chain(subnormal)
    {
        xs.push(values::midpoint_scan(x0, eval, 1 << 12));
    }

    // Close arguments can give the same midpoint
    xs.sort_by(f64::total_cmp);
    xs.dedup();
    xs
}
