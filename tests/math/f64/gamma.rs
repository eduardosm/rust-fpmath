use super::{MIN_MAX_ERROR, RUG_PREC, check_result, values};

#[test]
fn test_gamma() {
    let mut max_error: f64 = 0.0;
    test_with(Func::Gamma, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).gamma();
        let actual = fpmath::gamma(x);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_ln_gamma() {
    let mut max_error: f64 = 0.0;
    test_with(Func::LnGamma, |x| {
        let (expected, expected_sign) = if x.is_nan() || x == f64::NEG_INFINITY {
            (rug::Float::with_val(RUG_PREC, rug::float::Special::Nan), 0)
        } else {
            let (expected, ord) = rug::Float::with_val(RUG_PREC, x).ln_abs_gamma();
            let expected_sign = if x < 0.0 && x.fract() == 0.0 {
                0
            } else {
                ord as i8
            };
            (expected, expected_sign)
        };
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

fn test_with(func: Func, mut f: impl FnMut(f64)) {
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

    // Test a set of normal values at each binade
    for x in values::binades(-1022..=1023, 300) {
        f(x);
        f(-x);
    }

    // Test values evenly spread over the range where `gamma` is finite and not
    // zero
    for x in values::interval(-186.0, 172.0, 1_000_000) {
        f(x);
    }

    // Test more values in the problematic range
    for x in values::binades(1..=2, 250_000) {
        f(-x);
    }

    // Test the arguments around the negative integers (poles), whose results
    // are large or overflow
    for n in 1..=200 {
        for x in values::around(-f64::from(n), 100) {
            f(x);
        }
    }

    // Test the arguments around the boundaries of the argument reduction
    // (`n + 0.5`)
    for n in -186..=172 {
        for x in values::around(f64::from(n) + 0.5, 2) {
            f(x);
        }
    }

    // Test the arguments around the limits of the range
    for x in limits() {
        for x in values::around(x, 10_000) {
            f(x);
        }
    }

    // Test the arguments around the local extrema of `gamma`
    for x in &extrema {
        for x in values::around(x.to_f64(), 200) {
            f(x);
        }
    }

    // Test the arguments around the zeros of `ln_gamma` (where `|Γ(x)| = 1`),
    // whose results with `ln_gamma` are close to zero
    let zeros = ln_gamma_zeros(&extrema);
    let n = if func == Func::LnGamma { 10_000 } else { 1000 };
    for &x in &zeros {
        for x in values::around(x, n) {
            f(x);
        }
    }

    // Test the arguments around the zeros of `ln_gamma` with negative
    // arguments where its results are `±2^-16`, `±2^-7` and `±2^-4`, the
    // limits where it changes the evaluation method
    if func == Func::LnGamma {
        for x in ln_gamma_zero_limits(&zeros) {
            for x in values::around(x, 100) {
                f(x);
            }
        }
    }

    // Test the arguments further away from the zeros of `ln_gamma` with
    // negative arguments, from 2^14 to 2^48 ULPs, where it changes the
    // evaluation method at different distances
    if func == Func::LnGamma {
        for &x0 in zeros.iter().filter(|x| **x < 0.0) {
            for shift in 14..48 {
                for v in values::spread(16, shift) {
                    let bump = (1 << shift) | (v & ((1 << shift) - 1));
                    f(f64::from_bits(x0.to_bits() + bump));
                    f(f64::from_bits(x0.to_bits() - bump));
                }
            }
        }

        // Hard-to-round cases close to the zeros
        let hard_cases = [
            0xC003A7DC89BBC04B, // -2.456963611643038
            0xC003A7FBE7319E44, // -2.4570234357850342
            0xC003A7FC91A9B4CA, // -2.457024705878756
            0xC003A7FD93C0E770, // -2.457026628803696
            0xC005FB410BD22679, // -2.7476826594836976
            0xC005FB410C2B1733, // -2.747682662072202
            0xC005FB41122A6FCC, // -2.7476827067566543
            0xC005FB4318131FD9, // -2.747686565478926
        ];
        for bits in hard_cases {
            f(f64::from_bits(bits));
        }
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round
    for x in near_midpoints(func, &extrema) {
        f(x);
    }

    // Test integers
    for i in 1..=500 {
        let x = i as f64;
        f(x);
        f(-x);
    }
}

/// Returns the arguments at the limits of the range: where the functions
/// change the evaluation method (the Stirling series from 4 with `ln_gamma`
/// and from 6 with `gamma`, the reflection formula below -4 with `ln_gamma`
/// and from -10 with `gamma`, the approximations for tiny arguments below
/// 2^-60 with `ln_gamma` and 2^-64 with `gamma`, and the simpler Stirling
/// series from 2^52 with `ln_gamma`), where `gamma` overflows (its result is
/// the midpoint between the largest finite value and infinity) with large and
/// tiny arguments, and where `ln_gamma` overflows.
fn limits() -> [f64; 13] {
    let overflow: rug::Float =
        (rug::Float::with_val(EXT_PREC, 1) << 1024) - (rug::Float::with_val(EXT_PREC, 1) << 970);
    // Γ(x) ~= 1/x with tiny x
    let tiny = rug::Float::with_val(EXT_PREC, overflow.recip_ref()).to_f64();
    [
        4.0,
        6.0,
        -4.0,
        -10.0,
        fpmath::scalbn(1.0, -60),
        -fpmath::scalbn(1.0, -60),
        fpmath::scalbn(1.0, -64),
        -fpmath::scalbn(1.0, -64),
        fpmath::scalbn(1.0, 52),
        Func::Gamma.solve(&overflow, 171.6).to_f64(),
        tiny,
        -tiny,
        Func::LnGamma.solve(&overflow, 2.5e305).to_f64(),
    ]
}

/// Returns the arguments of the local extrema of `gamma` (where `ψ = Γ'/Γ` is
/// zero): its minimum with positive arguments (first), and one between each
/// pair of consecutive negative integers, down to -186.
fn gamma_extrema() -> Vec<rug::Float> {
    // ψ increases from -inf to +inf between consecutive poles, so its zeros can
    // be found by bisection (precise enough with 80 steps)
    let bisect = |lo: i32, hi: i32| {
        const PREC: u32 = 128;
        let mut a = rug::Float::with_val(PREC, lo) + (rug::Float::with_val(PREC, 1) >> 64);
        let mut b = rug::Float::with_val(PREC, hi) - (rug::Float::with_val(PREC, 1) >> 64);
        for _ in 0..80 {
            let m = rug::Float::with_val(PREC, &a + &b) / 2;
            if rug::Float::with_val(PREC, &m).digamma().is_sign_negative() {
                a = m;
            } else {
                b = m;
            }
        }
        a
    };
    [bisect(1, 2)]
        .into_iter()
        .chain((0..186).map(|n| bisect(-n - 1, -n)))
        .collect()
}

/// Returns the zeros of `ln_gamma` (where `|Γ(x)| = 1`): 1, 2, and two
/// between each pair of consecutive negative integers from -2 (on both sides
/// of the minimum of `|Γ|`, which is below one), until they are rounded to the
/// integers.
fn ln_gamma_zeros(extrema: &[rug::Float]) -> Vec<f64> {
    // ln(|Γ|) is positive close to the poles and negative at the minimum of
    // `|Γ|`, so its zeros can be found by bisection
    let bisect = |mut pos: rug::Float, mut neg: rug::Float| {
        for _ in 0..200 {
            let m = rug::Float::with_val(EXT_PREC, &pos + &neg) / 2;
            if Func::LnGamma.eval(&m).is_sign_negative() {
                neg = m;
            } else {
                pos = m;
            }
        }
        pos.to_f64()
    };

    let mut zeros = vec![1.0, 2.0];
    for x_e in extrema.iter().skip(1) {
        let n = x_e.to_f64().ceil();
        if n > -2.0 {
            continue;
        }
        let tiny = rug::Float::with_val(EXT_PREC, 1) >> 200;
        let left = bisect(rug::Float::with_val(EXT_PREC, n - 1.0) + &tiny, x_e.clone());
        let right = bisect(rug::Float::with_val(EXT_PREC, n) - tiny, x_e.clone());
        let len = zeros.len();
        zeros.extend([left, right].into_iter().filter(|x| x.fract() != 0.0));
        if zeros.len() == len {
            break;
        }
    }
    zeros
}

/// Returns the arguments around the zeros of `ln_gamma` with negative
/// arguments where its results are `±2^-16`, `±2^-7` and `±2^-4`
/// (calculated with Newton's method from the zeros).
fn ln_gamma_zero_limits(zeros: &[f64]) -> Vec<f64> {
    let mut xs = Vec::new();
    for &x0 in zeros.iter().filter(|x| **x < 0.0) {
        for e in [-16, -7, -4] {
            for sign in [1.0, -1.0] {
                let y = rug::Float::with_val(EXT_PREC, sign * fpmath::scalbn(1.0, e));
                let x = Func::LnGamma.solve(&y, x0).to_f64();
                // Newton's method could have converged to a solution between
                // other negative integers
                if x.is_finite() && x.floor() == x0.floor() {
                    xs.push(x);
                }
            }
        }
    }
    xs
}

/// Returns arguments whose results are very close to a midpoint between two
/// consecutive values.
///
/// Around the local extrema of `gamma` (which are also extrema of
/// `ln_gamma`), the functions are flat (their results change much less than
/// an ULP between consecutive arguments), so the arguments closest to the
/// inverse of a midpoint (calculated with Newton's method from a close
/// argument) have such results (see `values::midpoint_inverse`).
fn near_midpoints(func: Func, extrema: &[rug::Float]) -> Vec<f64> {
    let eval = |x: &rug::Float| func.eval(x);

    let mut xs = Vec::new();
    // The functions are flat closer to the extrema as their magnitude grows,
    // so the arguments are closer to the negative extrema, and there are more
    // of them around the minimum with positive arguments
    for (i, x_e) in extrema.iter().enumerate().filter(|(_, x)| **x > -40) {
        let (max_e, n) = if i == 0 { (-3, 100) } else { (-10, 4) };
        let dists = (-20..=max_e)
            .flat_map(|e| values::interval(fpmath::scalbn(1.0, e), fpmath::scalbn(1.0, e + 1), n));
        for d in dists {
            for x0 in [x_e.clone() - d, x_e.clone() + d] {
                let x0 = x0.to_f64();
                let x = values::midpoint_inverse(x0, eval, |y| func.solve(y, x0));
                // Newton's method could have converged to another solution
                if (x - x0).abs() <= d {
                    xs.extend(values::around(x, 1));
                }
            }
        }
    }

    // Close arguments can give the same midpoint
    xs.sort_by(f64::total_cmp);
    xs.dedup();
    xs
}
