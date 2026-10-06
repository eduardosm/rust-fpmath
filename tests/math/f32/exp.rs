use super::{MIN_MAX_ERROR, check_result, values};

#[test]
fn test_exp() {
    let mut max_error: f32 = 0.0;
    test_with(Func::Exp, |x| {
        let expected = fpmath::exp(f64::from(x));
        let actual = fpmath::exp(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_exp_m1() {
    let mut max_error: f32 = 0.0;
    test_with(Func::ExpM1, |x| {
        let expected = fpmath::exp_m1(f64::from(x));
        let actual = fpmath::exp_m1(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_exp2() {
    let mut max_error: f32 = 0.0;
    test_with(Func::Exp2, |x| {
        let expected = fpmath::exp2(f64::from(x));
        let actual = fpmath::exp2(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_exp10() {
    let mut max_error: f32 = 0.0;
    test_with(Func::Exp10, |x| {
        let expected = fpmath::exp10(f64::from(x));
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

impl Func {
    /// Calculates the function in `f64`.
    fn eval(self, x: f64) -> f64 {
        match self {
            Self::Exp => fpmath::exp(x),
            Self::ExpM1 => fpmath::exp_m1(x),
            Self::Exp2 => fpmath::exp2(x),
            Self::Exp10 => fpmath::exp10(x),
        }
    }

    /// Calculates the inverse of the function in `f64`.
    fn inverse(self, y: f64) -> f64 {
        match self {
            Self::Exp => fpmath::ln(y),
            Self::ExpM1 => fpmath::ln_1p(y),
            Self::Exp2 => fpmath::log2(y),
            Self::Exp10 => fpmath::log10(y),
        }
    }

    /// Returns `log_b(2)`, where `b` is the base of the exponential.
    fn log_b_2(self) -> f64 {
        match self {
            Self::Exp | Self::ExpM1 => std::f64::consts::LN_2,
            Self::Exp2 => 1.0,
            Self::Exp10 => std::f64::consts::LOG10_2,
        }
    }

    /// Returns the range of arguments with non-trivial results. Below it,
    /// the results are rounded to zero (or to -1 for `exp_m1`), and above it,
    /// they overflow.
    fn range(self) -> (f32, f32) {
        let lo = match self {
            Self::ExpM1 => exp2i(-30) - 1.0,
            _ => exp2i(-151),
        };
        (self.inverse(lo) as f32, self.inverse(exp2i(129)) as f32)
    }

    /// Returns the range of arguments with subnormal results, if any.
    fn subnormal_range(self) -> Option<(f32, f32)> {
        match self {
            Self::ExpM1 => None,
            _ => Some((
                self.inverse(exp2i(-150)) as f32,
                self.inverse(exp2i(-126)) as f32,
            )),
        }
    }

    /// Returns the results at the limits of the range: the midpoint between
    /// the largest finite value and infinity (above which the results
    /// overflow) and, at the other end, the smallest normal value, the
    /// smallest subnormal value and half of it (below which the results are
    /// rounded to zero), or the midpoint between -1 and the next value (below
    /// which the results of `exp_m1` are rounded to -1).
    fn limit_results(self) -> Vec<f64> {
        let overflow = exp2i(128) - exp2i(103);
        match self {
            Self::ExpM1 => vec![overflow, exp2i(-25) - 1.0],
            _ => vec![overflow, exp2i(-126), exp2i(-149), exp2i(-150)],
        }
    }
}

fn test_with(func: Func, mut f: impl FnMut(f32)) {
    // Test special values
    for x in values::specials() {
        f(x);
        f(-x);
    }

    // Test a set of subnormal values, whose results are trivial (one, or `x`
    // for `exp_m1`)
    for x in values::subnormals(10_000) {
        f(x);
        f(-x);
    }

    // Test a set of normal values at each binade, most of them have trivial
    // results (zero, one, infinity, `x` or -1)
    for x in values::binades(-126..=127, 1000) {
        f(x);
        f(-x);
    }

    // Test more values with small magnitudes, whose results are close to one
    // (or to `x` for `exp_m1`)
    for x in values::binades(-30..=-1, 200_000) {
        f(x);
        f(-x);
    }

    // Test values evenly spread over the whole range of non-trivial results,
    // so all the result binades are tested evenly
    let (lo, hi) = func.range();
    for x in values::interval(lo, hi, 20_000_000) {
        f(x);
    }

    // Exhaustive test of all the arguments with subnormal results
    if let Some((lo, hi)) = func.subnormal_range() {
        for x in values::interval_full(lo, hi) {
            f(x);
        }
    }

    // Test the arguments around the limits of the range
    for y in func.limit_results() {
        for x in values::around(func.inverse(y) as f32, 100_000) {
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
    for i in 1..=160 {
        let x = i as f32;
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
fn reduction_boundaries(func: Func) -> impl Iterator<Item = f32> {
    let c = func.log_b_2();
    let n = f64::from(REDUCTION_N);
    let (lo, hi) = func.range();
    let m_lo = (f64::from(lo) / c * n).floor() as i32;
    let m_hi = (f64::from(hi) / c * n).ceil() as i32;
    (m_lo..=m_hi).flat_map(move |m| {
        let power_of_two = (m % REDUCTION_N == 0).then(|| (f64::from(m / REDUCTION_N) * c) as f32);
        let boundary = (f64::from(2 * m + 1) * c / (2.0 * n)) as f32;
        power_of_two.into_iter().chain([boundary])
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
fn near_midpoints(func: Func) -> Vec<f32> {
    let small = values::binades(-27..=-2, 5000).flat_map(|x| [x, -x]);

    let eval = |x: f64| func.eval(x);
    let inverse = |y: f64| func.inverse(y);

    let mut xs = Vec::new();
    match func {
        Func::ExpM1 => {
            for x0 in small {
                // g(x) = exp_m1(x) - x, g'(x) = exp_m1(x)
                if let Some(x) = values::offset_midpoint(x0, fpmath::exp_m1, fpmath::exp_m1) {
                    xs.extend(values::around(x, 1));
                }
            }
            for x0 in values::interval(-17.0, -6.0, 100_000) {
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
    xs.sort_by(f32::total_cmp);
    xs.dedup();
    xs
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
