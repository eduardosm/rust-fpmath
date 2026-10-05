use super::{MIN_MAX_ERROR, RUG_PREC, check_result, values};

#[test]
fn test_exp() {
    let mut max_error: f64 = 0.0;
    test_with(Func::Exp, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).exp();
        let actual = fpmath::exp(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_exp_m1() {
    let mut max_error: f64 = 0.0;
    test_with(Func::ExpM1, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).exp_m1();
        let actual = fpmath::exp_m1(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_exp2() {
    let mut max_error: f64 = 0.0;
    test_with(Func::Exp2, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).exp2();
        let actual = fpmath::exp2(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_exp10() {
    let mut max_error: f64 = 0.0;
    test_with(Func::Exp10, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).exp10();
        let actual = fpmath::exp10(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

/// A function of the `exp` family.
#[derive(Copy, Clone)]
enum Func {
    Exp,
    ExpM1,
    Exp2,
    Exp10,
}

/// Precision used to calculate the test values.
const EXT_PREC: u32 = 256;

impl Func {
    /// Calculates the function with high precision.
    fn eval(self, x: &rug::Float) -> rug::Float {
        let x = rug::Float::with_val(EXT_PREC, x);
        match self {
            Self::Exp => x.exp(),
            Self::ExpM1 => x.exp_m1(),
            Self::Exp2 => x.exp2(),
            Self::Exp10 => x.exp10(),
        }
    }

    /// Calculates the inverse of the function with high precision.
    fn inverse(self, y: &rug::Float) -> rug::Float {
        let y = rug::Float::with_val(EXT_PREC, y);
        match self {
            Self::Exp => y.ln(),
            Self::ExpM1 => y.ln_1p(),
            Self::Exp2 => y.log2(),
            Self::Exp10 => y.log10(),
        }
    }

    /// Returns `log_b(2)`, where `b` is the base of the exponential. The
    /// arguments are reduced to `k * log_b(2) + r`, with an integer `k` and
    /// `|r| <= log_b(2) / 2`.
    fn log_b_2(self) -> rug::Float {
        let two = rug::Float::with_val(EXT_PREC, 2);
        match self {
            Self::Exp | Self::ExpM1 => two.ln(),
            Self::Exp2 => rug::Float::with_val(EXT_PREC, 1),
            Self::Exp10 => two.log10(),
        }
    }

    /// Returns the range of arguments with non-trivial results. Below it,
    /// the results are rounded to zero (or to -1 for `exp_m1`), and above it,
    /// they overflow.
    fn range(self) -> (f64, f64) {
        let lo = match self {
            Self::ExpM1 => exp2i(-60) - 1,
            _ => exp2i(-1076),
        };
        (
            self.inverse(&lo).to_f64(),
            self.inverse(&exp2i(1025)).to_f64(),
        )
    }

    /// Returns the range of arguments with subnormal results, if any.
    fn subnormal_range(self) -> Option<(f64, f64)> {
        match self {
            Self::ExpM1 => None,
            _ => Some((
                self.inverse(&exp2i(-1075)).to_f64(),
                self.inverse(&exp2i(-1022)).to_f64(),
            )),
        }
    }

    /// Returns the results at the limits of the range: the midpoint between
    /// the largest finite value and infinity (above which the results
    /// overflow) and, at the other end, the smallest normal value, the
    /// smallest subnormal value and half of it (below which the results are
    /// rounded to zero), or the midpoint between -1 and the next value (below
    /// which the results of `exp_m1` are rounded to -1).
    fn limit_results(self) -> Vec<rug::Float> {
        let overflow = exp2i(1024) - exp2i(970);
        match self {
            Self::ExpM1 => vec![overflow, exp2i(-54) - 1],
            _ => vec![overflow, exp2i(-1022), exp2i(-1074), exp2i(-1075)],
        }
    }
}

fn test_with(func: Func, mut f: impl FnMut(f64)) {
    // Test special values
    for x in values::specials() {
        f(x);
        f(-x);
    }

    // Test a set of subnormal values, whose results are trivial (one, or `x`
    // for `exp_m1`)
    for x in values::subnormals(1000) {
        f(x);
        f(-x);
    }

    // Test a set of normal values at each binade, most of them have trivial
    // results (zero, one, infinity, `x` or -1)
    for x in values::binades(-1022..=1023, 100) {
        f(x);
        f(-x);
    }

    // Test more values with small magnitudes, whose results are close to one
    // (or to `x` for `exp_m1`)
    for x in values::binades(-60..=-1, 20_000) {
        f(x);
        f(-x);
    }

    // Test values evenly spread over the whole range of non-trivial results,
    // so all the result binades are tested evenly
    let (lo, hi) = func.range();
    for x in values::interval(lo, hi, 2_000_000) {
        f(x);
    }

    // Test more values with subnormal results
    if let Some((lo, hi)) = func.subnormal_range() {
        for x in values::interval(lo, hi, 500_000) {
            f(x);
        }
    }

    // Test the arguments around the limits of the range
    for y in func.limit_results() {
        for x in values::around(func.inverse(&y).to_f64(), 10_000) {
            f(x);
        }
    }

    // Test the arguments around the boundaries of the argument reduction
    for x in reduction_boundaries(func) {
        for x in values::around(x, 2) {
            f(x);
        }
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round
    for x in near_midpoints(func) {
        f(x);
    }

    // Test integers
    for i in 1..=1100 {
        let x = i as f64;
        f(x);
        f(-x);
    }
}

/// Returns the arguments `k * log_b(2)` (with results close to powers of two)
/// and `(k + 1/2) * log_b(2)` (where `k` changes in the argument reduction)
/// for all `k` in the range of non-trivial results.
fn reduction_boundaries(func: Func) -> impl Iterator<Item = f64> {
    let c = func.log_b_2();
    let (lo, hi) = func.range();
    let k_lo = (lo / c.to_f64()).floor() as i32;
    let k_hi = (hi / c.to_f64()).ceil() as i32;
    (2 * k_lo..=2 * k_hi).map(move |k2| {
        let x: rug::Float = rug::Float::with_val(EXT_PREC, k2) * &c / 2;
        x.to_f64()
    })
}

/// Returns arguments whose results are very close to a midpoint between two
/// consecutive values.
///
/// Where the function is flat (its result changes much less than an ULP
/// between consecutive arguments), the arguments closest to the inverse of a
/// midpoint have such results. This happens with small arguments (results
/// close to one) and, for `exp_m1`, with large negative arguments (results
/// close to -1).
///
/// With small arguments, `exp_m1(x) = x + g(x)` is not flat, but `g` is. When
/// the result is in the same binade as `x`, the arguments closest to the
/// inverse of `g` at an odd multiple of half an ULP of `x` have results close
/// to a midpoint.
fn near_midpoints(func: Func) -> Vec<f64> {
    let small = values::binades(-56..=-2, 1000).flat_map(|x| [x, -x]);

    let eval = |x: &rug::Float| func.eval(x);
    let inverse = |y: &rug::Float| func.inverse(y);

    let mut xs = Vec::new();
    match func {
        Func::ExpM1 => {
            for x0 in small {
                // g(x) = exp_m1(x) - x, g'(x) = exp_m1(x)
                if let Some(x) = values::offset_midpoint(x0, eval, eval) {
                    xs.extend(values::around(x, 1));
                }
            }
            for x0 in values::interval(-37.0, -8.0, 100_000) {
                xs.extend(values::around(
                    values::midpoint_inverse(x0, eval, inverse),
                    1,
                ));
            }
        }
        _ => {
            for x0 in small {
                xs.extend(values::around(
                    values::midpoint_inverse(x0, eval, inverse),
                    1,
                ));
            }
        }
    }

    // Close arguments can give the same midpoint
    xs.sort_by(f64::total_cmp);
    xs.dedup();
    xs
}

/// Returns `2^e`.
fn exp2i(e: i32) -> rug::Float {
    rug::Float::with_val(EXT_PREC, 1) << e
}
