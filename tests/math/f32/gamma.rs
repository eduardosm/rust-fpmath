use super::{MIN_MAX_ERROR, check_result, values};

#[test]
fn test_gamma() {
    let mut max_error: f32 = 0.0;
    test_with(Func::Gamma, |x| {
        let expected = fpmath::gamma(f64::from(x));
        let actual = fpmath::gamma(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_ln_gamma() {
    let mut max_error: f32 = 0.0;
    test_with(Func::LnGamma, |x| {
        let (expected, expected_sign) = fpmath::ln_gamma(f64::from(x));
        let (actual, actual_sign) = fpmath::ln_gamma(x);

        assert_eq!(expected_sign, actual_sign);
        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

/// A function of the gamma family.
#[derive(Copy, Clone, PartialEq, Eq)]
enum Func {
    Gamma,
    LnGamma,
}

/// Precision used to calculate the test values.
const EXT_PREC: u32 = 256;

impl Func {
    /// Calculates the function (the logarithm of the absolute value of the
    /// gamma function for `ln_gamma`) in `f64`.
    fn eval_f64(self, x: f64) -> f64 {
        match self {
            Self::Gamma => fpmath::gamma(x),
            Self::LnGamma => fpmath::ln_gamma(x).0,
        }
    }

    /// Calculates the function (the logarithm of the absolute value of the
    /// gamma function for `ln_gamma`) with high precision.
    fn eval(self, x: &rug::Float) -> rug::Float {
        let x = rug::Float::with_val(EXT_PREC, x);
        match self {
            Self::Gamma => x.gamma(),
            Self::LnGamma => x.ln_abs_gamma().0,
        }
    }

    /// Calculates the derivative of the function with high precision, with
    /// `Γ'(x) = Γ(x) * ψ(x)` and `ln(|Γ(x)|)' = ψ(x)`.
    fn derivative(self, x: &rug::Float) -> rug::Float {
        let psi = rug::Float::with_val(EXT_PREC, x).digamma();
        match self {
            Self::Gamma => psi * self.eval(x),
            Self::LnGamma => psi,
        }
    }

    /// Solves `func(x) = y` with Newton's method, starting from `x0`.
    fn solve(self, y: f64, x0: f64) -> f64 {
        let mut x = rug::Float::with_val(EXT_PREC, x0);
        for _ in 0..100 {
            let dx: rug::Float = (self.eval(&x) - y) / self.derivative(&x);
            x -= &dx;
            if dx.is_zero() || dx.abs() < rug::Float::with_val(EXT_PREC, x.abs_ref()) >> 80 {
                break;
            }
        }
        x.to_f64()
    }
}

fn test_with(func: Func, mut f: impl FnMut(f32)) {
    let extrema = gamma_extrema();

    // Test special values
    for x in values::specials() {
        f(x);
        f(-x);
    }

    // Test a set of subnormal values
    for x in values::subnormals(100_000) {
        f(x);
        f(-x);
    }

    // Test across a wide range of normal numbers
    for x in values::binades(-126..=127, 5000) {
        f(x);
        f(-x);
    }

    // Test values evenly spread over the range where `gamma` is finite and not
    // zero
    for x in values::interval(-44.0, 36.0, 3_000_000) {
        f(x);
    }

    // Exhaustive test of all mantissas at binades that include the
    // problematic range
    for x in values::binades_full(1..=2) {
        f(-x);
    }

    // Test the arguments around the negative integers (poles), whose results
    // are large or overflow
    for n in 1..=50 {
        for x in values::around(-n as f32, 10_000) {
            f(x);
        }
    }

    // Test the arguments around the boundaries of the argument reduction
    // (`n + 0.5`)
    for n in -44..=36 {
        for x in values::around(n as f32 + 0.5, 2) {
            f(x);
        }
    }

    // Test the arguments around the limits of the range
    for x in limits() {
        for x in values::around(x, 100_000) {
            f(x);
        }
    }

    // Test the arguments around the local extrema of `gamma`
    for &x in extrema.iter() {
        for x in values::around(x as f32, 10_000) {
            f(x);
        }
    }

    // Test the arguments around the zeros of `ln_gamma` (where `|Γ(x)| = 1`),
    // whose results with `ln_gamma` are close to zero
    let zeros = ln_gamma_zeros(&extrema);
    let n = if func == Func::LnGamma {
        100_000
    } else {
        10_000
    };
    for &x in &zeros {
        for x in values::around(x, n) {
            f(x);
        }
    }

    // Test the arguments around the zeros of `ln_gamma` with negative
    // arguments where its results are `±2^-12`, the limit where it changes
    // the evaluation method
    if func == Func::LnGamma {
        for x in ln_gamma_zero_limits(&zeros) {
            for x in values::around(x, 1000) {
                f(x);
            }
        }
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round
    for x in near_midpoints(func, &extrema) {
        f(x);
    }

    // Test integers
    for i in 1..=500 {
        let x = i as f32;
        f(x);
        f(-x);
    }
}

/// Returns the arguments at the limits of the range: where the functions
/// change the evaluation method (the Stirling series from 4 with `ln_gamma`
/// and from 10 with `gamma`, and the reflection formula from -10), where
/// `gamma` overflows (its result is the midpoint between the largest finite
/// value and infinity) with large and tiny arguments, and where `ln_gamma`
/// overflows.
fn limits() -> [f32; 7] {
    let overflow = fpmath::scalbn(1.0, 128) - fpmath::scalbn(1.0, 103);
    // Γ(x) ~= 1/x with tiny x
    let tiny = (1.0 / overflow) as f32;
    [
        4.0,
        10.0,
        -10.0,
        Func::Gamma.solve(overflow, 35.0) as f32,
        tiny,
        -tiny,
        Func::LnGamma.solve(overflow, 4e36) as f32,
    ]
}

/// Returns the arguments around the zeros of `ln_gamma` with negative
/// arguments where its results are `±2^-12` (calculated with Newton's method
/// from the zeros).
fn ln_gamma_zero_limits(zeros: &[f32]) -> Vec<f32> {
    let mut xs = Vec::new();
    for &x0 in zeros.iter().filter(|x| **x < 0.0) {
        for y in [fpmath::scalbn(1.0, -12), -fpmath::scalbn(1.0, -12)] {
            let x = Func::LnGamma.solve(y, x0.into()) as f32;
            // Newton's method could have converged to a solution between
            // other negative integers
            if x.is_finite() && x.floor() == x0.floor() {
                xs.push(x);
            }
        }
    }
    xs
}

/// Returns the arguments of the local extrema of `gamma` (where `ψ = Γ'/Γ` is
/// zero): its minimum with positive arguments (first), and one between each
/// pair of consecutive negative integers, down to -50.
fn gamma_extrema() -> Vec<f64> {
    // ψ increases from -inf to +inf between consecutive poles, so its zeros can
    // be found by bisection (precise enough with 80 steps)
    let bisect = |lo: i32, hi: i32| {
        const PREC: u32 = 128;
        let mut a: rug::Float =
            rug::Float::with_val(PREC, lo) + (rug::Float::with_val(PREC, 1) >> 64);
        let mut b = rug::Float::with_val(PREC, hi) - (rug::Float::with_val(PREC, 1) >> 64);
        for _ in 0..80 {
            let m = rug::Float::with_val(PREC, &a + &b) / 2;
            if rug::Float::with_val(PREC, &m).digamma().is_sign_negative() {
                a = m;
            } else {
                b = m;
            }
        }
        a.to_f64()
    };
    [bisect(1, 2)]
        .into_iter()
        .chain((0..50).map(|n| bisect(-n - 1, -n)))
        .collect()
}

/// Returns the zeros of `ln_gamma` (where `|Γ(x)| = 1`): 1, 2, and two
/// between each pair of consecutive negative integers from -2 (on both sides
/// of the minimum of `|Γ|`, which is below one), until they are rounded to the
/// integers.
fn ln_gamma_zeros(extrema: &[f64]) -> Vec<f32> {
    // ln(|Γ|) is positive close to the poles and negative at the minimum of
    // `|Γ|`, so its zeros can be found by bisection
    let bisect = |mut pos: rug::Float, mut neg: rug::Float| {
        for _ in 0..100 {
            let m = rug::Float::with_val(EXT_PREC, &pos + &neg) / 2;
            if Func::LnGamma.eval(&m).is_sign_negative() {
                neg = m;
            } else {
                pos = m;
            }
        }
        pos.to_f32()
    };

    let mut zeros = vec![1.0, 2.0];
    for &x_e in extrema.iter().skip(1) {
        let n = x_e.ceil();
        if n > -2.0 {
            continue;
        }
        let tiny = rug::Float::with_val(EXT_PREC, 1) >> 200;
        let x_e = rug::Float::with_val(EXT_PREC, x_e);
        let left = bisect(rug::Float::with_val(EXT_PREC, n - 1.0) + &tiny, x_e.clone());
        let right = bisect(rug::Float::with_val(EXT_PREC, n) - tiny, x_e);
        let len = zeros.len();
        zeros.extend([left, right].into_iter().filter(|x| x.fract() != 0.0));
        if zeros.len() == len {
            break;
        }
    }
    zeros
}

/// Returns arguments whose results are very close to a midpoint between two
/// consecutive values.
///
/// Around the local extrema of `gamma` (which are also extrema of
/// `ln_gamma`), the functions are flat (their results change much less than
/// an ULP between consecutive arguments), so the arguments closest to the
/// inverse of a midpoint (calculated with Newton's method from a close
/// argument) have such results (see `values::midpoint_inverse`).
fn near_midpoints(func: Func, extrema: &[f64]) -> Vec<f32> {
    let eval = |x: f64| func.eval_f64(x);

    let mut xs = Vec::new();
    // The functions are flat closer to the extrema as their magnitude grows,
    // so the arguments are closer to the negative extrema, and there are more
    // of them around the minimum with positive arguments
    for (i, &x_e) in extrema.iter().enumerate().filter(|(_, x)| **x > -35.0) {
        let (max_e, n) = if i == 0 { (-3, 200) } else { (-6, 4) };
        let dists = (-10..=max_e)
            .flat_map(|e| values::interval(fpmath::scalbn(1.0, e), fpmath::scalbn(1.0, e + 1), n));
        for d in dists {
            for x0 in [x_e - f64::from(d), x_e + f64::from(d)] {
                let x0 = x0 as f32;
                let x = values::midpoint_inverse(x0, eval, |y| func.solve(y, x0.into()));
                // Newton's method could have converged to another solution
                if (f64::from(x) - f64::from(x0)).abs() <= f64::from(d) {
                    xs.extend(values::around(x, 1));
                }
            }
        }
    }

    // Close arguments can give the same midpoint
    xs.sort_by(f32::total_cmp);
    xs.dedup();
    xs
}
