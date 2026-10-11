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
fn test_exp2_m1() {
    let mut max_error: f64 = 0.0;
    test_with(Func::Exp2M1, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).exp2_m1();
        let actual = fpmath::exp2_m1(x);

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

#[test]
fn test_exp10_m1() {
    let mut max_error: f64 = 0.0;
    test_with(Func::Exp10M1, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).exp10_m1();
        let actual = fpmath::exp10_m1(x);

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
    Exp2M1,
    Exp10,
    Exp10M1,
}

/// Precision used to calculate the test values.
const EXT_PREC: u32 = 256;

impl Func {
    /// Returns whether the function is `b^x - 1` (instead of `b^x`).
    fn is_m1(self) -> bool {
        matches!(self, Self::ExpM1 | Self::Exp2M1 | Self::Exp10M1)
    }

    /// Calculates the function with high precision.
    fn eval(self, x: &rug::Float) -> rug::Float {
        let x = rug::Float::with_val(EXT_PREC, x);
        match self {
            Self::Exp => x.exp(),
            Self::ExpM1 => x.exp_m1(),
            Self::Exp2 => x.exp2(),
            Self::Exp2M1 => x.exp2_m1(),
            Self::Exp10 => x.exp10(),
            Self::Exp10M1 => x.exp10_m1(),
        }
    }

    /// Calculates the inverse of the function with high precision.
    fn inverse(self, y: &rug::Float) -> rug::Float {
        let y = rug::Float::with_val(EXT_PREC, y);
        match self {
            Self::Exp => y.ln(),
            Self::ExpM1 => y.ln_1p(),
            Self::Exp2 => y.log2(),
            Self::Exp2M1 => y.log2_1p(),
            Self::Exp10 => y.log10(),
            Self::Exp10M1 => y.log10_1p(),
        }
    }

    /// Returns `log_b(2)`, where `b` is the base of the exponential.
    fn log_b_2(self) -> rug::Float {
        let two = rug::Float::with_val(EXT_PREC, 2);
        match self {
            Self::Exp | Self::ExpM1 => two.ln(),
            Self::Exp2 | Self::Exp2M1 => rug::Float::with_val(EXT_PREC, 1),
            Self::Exp10 | Self::Exp10M1 => two.log10(),
        }
    }

    /// Returns the range of arguments with non-trivial results. Below it,
    /// the results are rounded to zero (or to -1 for `b^x - 1`), and above
    /// it, they overflow.
    fn range(self) -> (f64, f64) {
        let lo = if self.is_m1() {
            exp2i(-60) - 1
        } else {
            exp2i(-1076)
        };
        (
            self.inverse(&lo).to_f64(),
            self.inverse(&exp2i(1025)).to_f64(),
        )
    }

    /// Returns the range of arguments with subnormal results that are not
    /// trivial, if any. The results of `b^x - 1` are subnormal around zero,
    /// where they are `x * ln(b)` rounded (which is trivial for `exp_m1`).
    fn subnormal_range(self) -> Option<(f64, f64)> {
        match self {
            Self::ExpM1 => None,
            Self::Exp2M1 | Self::Exp10M1 => Some((
                self.inverse(&-exp2i(-1022)).to_f64(),
                self.inverse(&exp2i(-1022)).to_f64(),
            )),
            _ => Some((
                self.inverse(&exp2i(-1075)).to_f64(),
                self.inverse(&exp2i(-1022)).to_f64(),
            )),
        }
    }

    /// Returns the range of arguments where the results of `b^x - 1` are
    /// close to -1 and the function is flat (see [`near_midpoints`]).
    fn flat_range(self) -> Option<(f64, f64)> {
        match self {
            Self::ExpM1 => Some((-37.0, -8.0)),
            Self::Exp2M1 => Some((-53.5, -11.5)),
            Self::Exp10M1 => Some((-16.0, -3.5)),
            _ => None,
        }
    }

    /// Returns the results at the limits of the range: the midpoint between
    /// the largest finite value and infinity (above which the results
    /// overflow) and, at the other end, the smallest normal value, the
    /// smallest subnormal value and half of it (below which the results are
    /// rounded to zero), or the midpoint between -1 and the next value (below
    /// which the results of `b^x - 1` are rounded to -1).
    fn limit_results(self) -> Vec<rug::Float> {
        let overflow = exp2i(1024) - exp2i(970);
        if self.is_m1() {
            vec![overflow, exp2i(-54) - 1]
        } else {
            vec![overflow, exp2i(-1022), exp2i(-1074), exp2i(-1075)]
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
    // for `exp_m1`) or `x * ln(b)` rounded (for `exp2_m1` and `exp10_m1`)
    for x in values::subnormals(1000) {
        f(x);
        f(-x);
    }

    // Test a set of normal values at each binade, most of them have trivial
    // results (zero, one, infinity, `x` or -1) or `x * ln(b)` rounded
    for x in values::binades(-1022..=1023, 100) {
        f(x);
        f(-x);
    }

    // Test more values with small magnitudes, whose results are close to one
    // (or to `x * ln(b)` for `b^x - 1`)
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

/// Number of parts in which the argument reduction splits each power of two:
/// the arguments are reduced to `m * log_b(2) / N + r`, with an integer `m`
/// and `|r| <= log_b(2) / (2 * N)`.
const REDUCTION_N: i32 = 128;

/// Returns the arguments `k * log_b(2)` (with results close to powers of two)
/// and `(m + 1/2) * log_b(2) / N` (where `m` changes in the argument reduction)
/// for all `k` and `m` in the range of non-trivial results.
fn reduction_boundaries(func: Func) -> impl Iterator<Item = f64> {
    let c = func.log_b_2();
    let (lo, hi) = func.range();
    let m_lo = (lo / c.to_f64() * f64::from(REDUCTION_N)).floor() as i32;
    let m_hi = (hi / c.to_f64() * f64::from(REDUCTION_N)).ceil() as i32;
    (m_lo..=m_hi).flat_map(move |m| {
        let power_of_two = (m % REDUCTION_N == 0).then(|| {
            let x: rug::Float = rug::Float::with_val(EXT_PREC, m / REDUCTION_N) * &c;
            x.to_f64()
        });
        let boundary: rug::Float =
            rug::Float::with_val(EXT_PREC, 2 * m + 1) * &c / (2 * REDUCTION_N);
        power_of_two.into_iter().chain([boundary.to_f64()])
    })
}

/// Returns arguments whose results are very close to a midpoint between two
/// consecutive values.
///
/// Where the function is flat (its result changes much less than an ULP
/// between consecutive arguments), the arguments closest to the inverse of a
/// midpoint have such results. This happens with small arguments (results
/// close to one) and, for `b^x - 1`, with large negative arguments (results
/// close to -1).
///
/// With small arguments, `exp_m1(x) = x + g(x)` is not flat, but `g` is. When
/// the result is in the same binade as `x`, the arguments closest to the
/// inverse of `g` at an odd multiple of half an ULP of `x` have results close
/// to a midpoint.
///
/// `exp2_m1` and `exp10_m1` are close to `x * ln(b)` with small arguments, so
/// the arguments with results close to a midpoint are searched among the ones
/// close to other values (see `values::midpoint_scan`), also with tiny
/// arguments, down to the smallest ones with normal results.
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
        }
        Func::Exp2M1 | Func::Exp10M1 => {
            // Fewer values than for the other functions, each one needs
            // several evaluations
            let x0s = values::binades(-56..=-2, 100)
                .chain(values::binades((-1021..=-57).step_by(16), 20));
            for x0 in x0s.flat_map(|x| [x, -x]) {
                xs.push(values::midpoint_scan(x0, eval, 1 << 12));
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
    if let Some((lo, hi)) = func.flat_range() {
        for x0 in values::interval(lo, hi, 100_000) {
            xs.extend(values::around(
                values::midpoint_inverse(x0, eval, inverse),
                1,
            ));
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
