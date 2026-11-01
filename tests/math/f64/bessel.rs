use rug::ops::Pow as _;

use super::{ERROR_LIMIT, MIN_MAX_ERROR, RUG_PREC, calc_error_ulp, check_result, values};

#[test]
fn test_j0() {
    let mut max_error: f64 = 0.0;
    test_with(Func::J0, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).j0();
        let actual = fpmath::j0(x);
        assert_result_eq!(fpmath::j0(-x), actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_j1() {
    let mut max_error: f64 = 0.0;
    test_with(Func::J1, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).j1();
        let actual = fpmath::j1(x);
        assert_result_eq!(fpmath::j1(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_y0() {
    let mut max_error: f64 = 0.0;
    test_with(Func::Y0, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).y0();
        let actual = fpmath::y0(x);
        if x > 0.0 {
            assert!(fpmath::y0(-x).is_nan());
        }

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_y1() {
    let mut max_error: f64 = 0.0;
    test_with(Func::Y1, |x| {
        let expected = rug::Float::with_val(RUG_PREC, x).y1();
        let actual = fpmath::y1(x);
        if x > 0.0 {
            assert!(fpmath::y1(-x).is_nan());
        }

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_jn() {
    let mut max_error: f64 = 0.0;
    for n in ORDERS {
        test_order_with(Kind::First, n, |x| {
            let actual = fpmath::jn(n, x);
            // jn(-n, x) = (-1)^n * jn(n, x) and jn(n, -x) = (-1)^n * jn(n, x)
            let odd = if n % 2 == 0 { actual } else { -actual };
            assert_result_eq!(fpmath::jn(-n, x), odd);
            assert_result_eq!(fpmath::jn(n, -x), odd);
            assert_result_eq!(fpmath::jn(-n, -x), actual);

            check_order_result(Kind::First, n, x, actual, &mut max_error);
        });
    }
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_yn() {
    let mut max_error: f64 = 0.0;
    for n in ORDERS {
        test_order_with(Kind::Second, n, |x| {
            let actual = fpmath::yn(n, x);
            // yn(-n, x) = (-1)^n * yn(n, x)
            let odd = if n % 2 == 0 { actual } else { -actual };
            assert_result_eq!(fpmath::yn(-n, x), odd);
            if x > 0.0 {
                assert!(fpmath::yn(n, -x).is_nan());
            }

            check_order_result(Kind::Second, n, x, actual, &mut max_error);
        });
    }
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_rug_orders() {
    // The recurrences that `eval_order` uses must match the functions of
    // `rug`, also close to their zeros, where the terms of the recurrences
    // cancel
    for n in [2, 3, 5, 12, 40] {
        let nf = f64::from(n);
        let xs = values::interval(0.01, 4.0 * nf + 40.0, 200)
            .chain(order_zeros(Kind::First, n, 4))
            .chain(order_zeros(Kind::Second, n, 4));
        for x in xs {
            let xr = rug::Float::with_val(53, x);
            for kind in [Kind::First, Kind::Second] {
                let exact = match kind {
                    Kind::First => rug::Float::with_val(EXT_PREC, &xr).jn(n),
                    Kind::Second => rug::Float::with_val(EXT_PREC, &xr).yn(n),
                };
                let error = (eval_order(kind, n, &xr, RUG_PREC + 8) - &exact) / exact;
                assert!(
                    error.abs() < rug::Float::with_val(53, 1) >> RUG_PREC,
                    "n = {n}, x = {x:e}"
                );
            }
        }
    }
}

#[test]
fn test_large_order_reference() {
    // `eval_large_order` must match `eval_order` with the orders where both
    // can be used, including the arguments close to the order, where it uses
    // its recurrences
    for n in [1024, 4096, 4099] {
        for x in large_order_args(n) {
            let xr = rug::Float::with_val(53, x);
            let first = eval_order(Kind::First, n as i32, &xr, EXT_PREC);
            let second = eval_order(Kind::Second, n as i32, &xr, EXT_PREC);
            let modulus = first.clone().hypot(&second);
            for (kind, exact) in [(Kind::First, first), (Kind::Second, second)] {
                // Relative to the result, or to the modulus beyond the order
                let scale = if x <= f64::from(n) {
                    exact.clone().abs()
                } else {
                    modulus.clone()
                };
                let error = (eval_large_order(kind, n, &xr) - exact).abs() / scale;
                assert!(
                    error < rug::Float::with_val(53, 1) >> 100u32,
                    "n = {n}, x = {x:e}, kind = {kind:?}, error = {:e}",
                    error.to_f64()
                );
            }
        }
    }
}

#[test]
fn test_large_orders() {
    // From the first orders with which `fpmath` uses the asymptotic
    // expansions (also tested with `eval_order`, but with fewer arguments
    // close to the order), to the largest ones (`i32::MIN` for 2^31)
    let mut max_error: f64 = 0.0;
    for n in [
        1024,
        1031,
        4096,
        4099,
        16_384,
        50_001,
        1_000_003,
        (1 << 31) - 1,
        1 << 31,
    ] {
        let order = if n == 1 << 31 { i32::MIN } else { n as i32 };
        for x in large_order_args(n) {
            let xr = rug::Float::with_val(53, x);
            let first = eval_large_order(Kind::First, n, &xr);
            let second = eval_large_order(Kind::Second, n, &xr);
            let eval = |kind| match kind {
                Kind::First => first.clone(),
                Kind::Second => second.clone(),
            };
            let actual = fpmath::jn(order, x);
            check_order_result_with(Kind::First, n, x, actual, eval, &mut max_error);
            let actual = fpmath::yn(order, x);
            check_order_result_with(Kind::Second, n, x, actual, eval, &mut max_error);
        }

        // Results very close to a midpoint between two consecutive values,
        // below the order (where the error limit is the one of the other
        // functions), from about e^-700 (`jn`) and e^700 (`yn`) to where
        // `eval_large_order` uses its recurrences, and from there to the
        // order: the Debye expansions of `fpmath` close to their limit, and
        // Olver's beyond it
        let x_lo = order_core_arg(n, 700.0);
        let [x_ref, _] = order_dist_args(n, DEBYE_REF_DIST);
        let nf = f64::from(n);
        for kind in [Kind::First, Kind::Second] {
            let eval = |x: &rug::Float| eval_large_order(kind, n, x);
            let xs = values::interval(x_lo, x_ref, 16).chain(values::interval(x_ref, nf, 16));
            for x0 in xs {
                let x = values::midpoint_scan(x0, eval, 1 << 12);
                let actual = match kind {
                    Kind::First => fpmath::jn(order, x),
                    Kind::Second => fpmath::yn(order, x),
                };
                let xr = rug::Float::with_val(53, x);
                let eval = |kind| eval_large_order(kind, n, &xr);
                check_order_result_with(kind, n, x, actual, eval, &mut max_error);
            }
        }
    }
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

/// Returns the arguments `x < n` and `x > n` such that `R^3 / n^2 = d`,
/// where `R = sqrt(|n^2 - x^2|)`.
fn order_dist_args(n: u32, d: f64) -> [f64; 2] {
    let nf = f64::from(n);
    let c = fpmath::cbrt(d * nf * nf);
    let r2 = c * c;
    [fpmath::sqrt(nf * nf - r2), fpmath::sqrt(nf * nf + r2)]
}

/// Returns the argument `x < n` such that `n * (atanh(w) - w) = e`, with
/// `w = sqrt(1 - (x / n)^2)`, where `jn(n, x)` is about `exp(-e)` and
/// `yn(n, x)` about `-exp(e)`.
fn order_core_arg(n: u32, e: f64) -> f64 {
    let nf = f64::from(n);
    let (mut w_lo, mut w_hi): (f64, f64) = (0.0, 1.0);
    for _ in 0..100 {
        let w = 0.5 * (w_lo + w_hi);
        if nf * (fpmath::atanh(w) - w) < e {
            w_lo = w;
        } else {
            w_hi = w;
        }
    }
    nf * fpmath::sqrt(1.0 - w_lo * w_lo)
}

/// Returns the arguments to test with the large order `n`: around the limits
/// of the asymptotic expansions of `fpmath` (`R^3 / n^2 = 256`, where
/// `R = sqrt(|n^2 - x^2|)`) and of `eval_large_order`, spread between them,
/// where the results underflow and overflow, and around the limits where the
/// evaluation method changes beyond the order.
fn large_order_args(n: u32) -> Vec<f64> {
    let nf = f64::from(n);
    let [lo, hi] = order_dist_args(n, 256.0);
    let [ref_lo, ref_hi] = order_dist_args(n, DEBYE_REF_DIST);

    let mut xs = Vec::new();
    // Around the limits, where the methods of `fpmath` and of
    // `eval_large_order` change
    for x0 in [lo, hi, ref_lo, ref_hi, nf] {
        xs.extend(values::around(x0, 2));
    }
    xs.extend(values::interval(lo, hi, 40));
    xs.extend(values::interval(ref_lo, lo, 10));
    xs.extend(values::interval(hi, ref_hi, 10));

    // Below the order, down to where `jn` underflows and `yn` overflows
    for e in [745.0, 710.0] {
        let x0 = order_core_arg(n, e);
        xs.extend(values::around(x0, 3));
        xs.extend(values::interval(x0, ref_lo, 10));
    }

    // Beyond the order: where the phase is calculated differently (R = n),
    // spread up to the limit from which the functions are calculated from
    // their modulus and phase (n^2 * 2^16), and huge arguments
    let s2 = std::f64::consts::SQRT_2 * nf;
    xs.extend(values::around(s2, 2));
    xs.extend([1.01, 1.1, 1.5, 2.0, 10.0, 1e3].map(|f| f * nf));
    xs.extend(values::around(fpmath::scalbn(nf * nf, 16), 2));
    xs.extend([fpmath::scalbn(1.0, 200), fpmath::scalbn(1.0, 600), f64::MAX]);

    xs.retain(|&x| x > 0.0 && x.is_finite());
    xs
}

/// Number of terms of the Debye expansions of `eval_large_order`, which are
/// more than the ones of `fpmath` (see `DEBYE_REF_DIST`).
const DEBYE_REF_TERMS: usize = 20;

/// `eval_large_order` uses the Debye expansions when `R^3 / n^2` is at least
/// this, where `R = sqrt(|n^2 - x^2|)`. With `DEBYE_REF_TERMS` terms, the
/// first one that is left out is less than 2^-114 then. It is larger than
/// the limit of `fpmath`, so the results of `fpmath` around that limit (where
/// it changes from the Debye expansions to Olver's) are compared with the
/// recurrences of `eval_large_order`.
const DEBYE_REF_DIST: f64 = 512.0;

/// Precision of `eval_large_order`, which is enough for the reduction of
/// the phase of the arguments up to 2^1024.
const DEBYE_REF_PREC: u32 = 1200;

/// Returns the coefficients of the polynomials `W[k](p)` (`W[k][j]` for
/// `p^j`), for `1 <= k <= DEBYE_REF_TERMS`, such that
/// `ln(sum(U[k](p) / n^k)) = sum(W[k](p) / n^k)`, where `U[k]` are the
/// polynomials of the Debye expansions (DLMF 10.41.10), with `U[0] = 1`.
fn debye_coeffs() -> &'static [Vec<rug::Float>] {
    static COEFFS: std::sync::OnceLock<Vec<Vec<rug::Float>>> = std::sync::OnceLock::new();
    COEFFS.get_or_init(|| {
        let prec = 1024;
        let zero = || rug::Float::new(prec);
        // U[k + 1](p) = p^2 * (1 - p^2) / 2 * U[k]'(p)
        //     + 1/8 * integral(0..p, (1 - 5 * t^2) * U[k](t) dt)
        let mut u = vec![vec![rug::Float::with_val(prec, 1)]];
        for k in 0..DEBYE_REF_TERMS {
            let mut next = vec![zero(); u[k].len() + 3];
            for (j, c) in u[k].iter().enumerate() {
                let jf = j as u32;
                let t = rug::Float::with_val(prec, c * jf) / 2u32;
                next[j + 1] += &t;
                next[j + 3] -= &t;
                next[j + 1] += rug::Float::with_val(prec, c / (8 * (jf + 1)));
                next[j + 3] -= rug::Float::with_val(prec, c * 5u32) / (8 * (jf + 3));
            }
            u.push(next);
        }
        // The logarithm of S = sum(U[k] * e^k), from S' = (ln S)' * S:
        // k * W[k] = k * U[k] - sum(j * W[j] * U[k - j]), for 1 <= j < k
        let mut w: Vec<Vec<rug::Float>> = vec![Vec::new()];
        for k in 1..=DEBYE_REF_TERMS {
            let mut acc: Vec<rug::Float> = u[k]
                .iter()
                .map(|c| rug::Float::with_val(prec, c * k as u32))
                .collect();
            for j in 1..k {
                for (ia, ca) in w[j].iter().enumerate() {
                    for (ib, cb) in u[k - j].iter().enumerate() {
                        acc[ia + ib] -= rug::Float::with_val(prec, ca * cb) * j as u32;
                    }
                }
            }
            w.push(acc.into_iter().map(|c| c / k as u32).collect());
        }
        w
    })
}

/// Calculates the Bessel function of the kind `kind` and order `n` of
/// `x > 0` from the Debye expansions, for the large orders with which
/// `eval_order` is too slow.
///
/// They are not valid close to `x = n`, where it uses the recurrences of
/// `eval_order` from the closest orders where they are, in the direction in
/// which they are stable: backward from orders above `x` with the functions
/// of the first kind, and forward from orders below it with the ones of the
/// second kind.
fn eval_large_order(kind: Kind, n: u32, x: &rug::Float) -> rug::Float {
    let xf = x.to_f64();
    let valid = |m: u32| {
        let m = f64::from(m);
        let r2 = ((m - xf) * (m + xf)).abs();
        r2 * r2 * r2 >= DEBYE_REF_DIST * DEBYE_REF_DIST * ((m * m) * (m * m))
    };
    if valid(n) {
        return debye_eval(kind, n, x);
    }

    let prec = DEBYE_REF_PREC;
    let xr = rug::Float::with_val(prec, x);
    // a[k + 1] = (2 * k / x) * a[k] - a[k - 1], in either direction
    let step = |a: &rug::Float, k: u32| rug::Float::with_val(prec, a * (2 * u64::from(k))) / &xr;
    match kind {
        Kind::First => {
            let mut m = n.max(xf as u32 + 1);
            while !valid(m) {
                m += 1;
            }
            let (mut a, mut b) = (debye_eval(kind, m + 1, x), debye_eval(kind, m, x));
            for k in (n + 1..=m).rev() {
                let c = step(&b, k) - &a;
                a = b;
                b = c;
            }
            b
        }
        Kind::Second => {
            let mut m = n.min(fpmath::ceil(xf) as u32 - 1);
            while !valid(m) {
                m -= 1;
            }
            let (mut a, mut b) = (debye_eval(kind, m - 1, x), debye_eval(kind, m, x));
            for k in m..n {
                let c = step(&b, k) - &a;
                a = b;
                b = c;
            }
            b
        }
    }
}

/// Calculates the Bessel function of the kind `kind` and order `n` of `x > 0`
/// from the Debye expansions (DLMF 10.19), where they are valid (see
/// `eval_large_order`).
///
/// With `R = sqrt(|n^2 - x^2|)`, `ln(sum(U[k](p) / n^k)) = sum(W[k](p) / n^k)`
/// (see `debye_coeffs`) and `S(p, n) = exp(sum(W[k](p) / n^k))`:
///
/// * `x < n`, with `w = R / n`:
///   * `J(x) = exp(-n * (atanh(w) - w)) / sqrt(2 * π * R) * S(1 / w, n)`
///   * `Y(x) = -exp(n * (atanh(w) - w)) / sqrt(π * R / 2) * S(1 / w, -n)`
/// * `x > n`, with `v = R / n`:
///   `J(x) + i * Y(x) = sqrt(2 / (π * R)) * exp(i * (n * (v - atan(v)) - π/4)) * S(i / v, -n)`
fn debye_eval(kind: Kind, n: u32, x: &rug::Float) -> rug::Float {
    let prec = DEBYE_REF_PREC;
    let pi = rug::Float::with_val(prec, rug::float::Constant::Pi);
    let nr = rug::Float::with_val(prec, n);
    let xr = rug::Float::with_val(prec, x);
    let n2 = rug::Float::with_val(prec, nr.square_ref());
    let x2 = rug::Float::with_val(prec, xr.square_ref());
    let below = xr < nr;
    let r = if below { n2 - x2 } else { x2 - n2 }.sqrt();
    let t = rug::Float::with_val(prec, &r / &nr);
    let q = rug::Float::with_val(prec, t.recip_ref());

    // sum(W[k](p) / m^k), with p = q (x < n) or p = i * q, as `re + i * im`
    let m = match (below, kind) {
        (true, Kind::First) => nr.clone(),
        _ => -nr.clone(),
    };
    let (mut re, mut im) = (rug::Float::new(prec), rug::Float::new(prec));
    for (k, w) in debye_coeffs().iter().enumerate().skip(1) {
        let scale = rug::Float::with_val(prec, (&m).pow(k as u32)).recip();
        for (j, c) in w.iter().enumerate() {
            let term = rug::Float::with_val(prec, (&q).pow(j as u32)) * c * &scale;
            // i^j when p = i * q
            match if below { 0 } else { j % 4 } {
                0 => re += term,
                1 => im += term,
                2 => re -= term,
                _ => im -= term,
            }
        }
    }

    if below {
        let core = (rug::Float::with_val(prec, t.atanh_ref()) - &t) * &nr;
        match kind {
            Kind::First => {
                let d = rug::Float::with_val(prec, &r * &pi) * 2u32;
                (re - core).exp() / d.sqrt()
            }
            Kind::Second => {
                let d = rug::Float::with_val(prec, &r * &pi) / 2u32;
                -(re + core).exp() / d.sqrt()
            }
        }
    } else {
        let theta = (rug::Float::with_val(prec, t.atan_ref()) - &t) * -nr - pi.clone() / 4u32 + im;
        let modulus = (rug::Float::with_val(prec, 2u32) / (r * pi)).sqrt() * re.exp();
        match kind {
            Kind::First => modulus * theta.cos(),
            Kind::Second => modulus * theta.sin(),
        }
    }
}

/// A Bessel function of the first (`J`) or second (`Y`) kind, of order 0 or 1.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Func {
    J0,
    J1,
    Y0,
    Y1,
}

/// Precision used to calculate the test values.
const EXT_PREC: u32 = 256;

impl Func {
    fn order(self) -> u32 {
        match self {
            Self::J0 | Self::Y0 => 0,
            Self::J1 | Self::Y1 => 1,
        }
    }

    /// Calculates the function with precision `prec`.
    fn eval_prec(self, x: &rug::Float, prec: u32) -> rug::Float {
        let x = rug::Float::with_val(prec, x);
        match self {
            Self::J0 => x.j0(),
            Self::J1 => x.j1(),
            Self::Y0 => x.y0(),
            Self::Y1 => x.y1(),
        }
    }

    /// Calculates the function with high precision.
    fn eval(self, x: &rug::Float) -> rug::Float {
        self.eval_prec(x, EXT_PREC)
    }

    /// Calculates the derivative of the function with precision `prec`, with
    /// `C0'(x) = -C1(x)` and `C1'(x) = C0(x) - C1(x) / x`.
    fn derivative_prec(self, x: &rug::Float, prec: u32) -> rug::Float {
        let x = rug::Float::with_val(prec, x);
        match self {
            Self::J0 => -x.j1(),
            Self::J1 => x.clone().j0() - x.clone().j1() / x,
            Self::Y0 => -x.y1(),
            Self::Y1 => x.clone().y0() - x.clone().y1() / x,
        }
    }

    /// Calculates the second derivative of the function with precision
    /// `prec`, with `C''(x) = -C'(x) / x - (1 - order^2 / x^2) * C(x)`.
    fn derivative2_prec(self, x: &rug::Float, prec: u32) -> rug::Float {
        let x = rug::Float::with_val(prec, x);
        let k = 1u32 - self.order().pow(2) / rug::Float::with_val(prec, x.square_ref());
        -self.derivative_prec(&x, prec) / &x - k * self.eval_prec(&x, prec)
    }

    /// Solves `func(x) = y` with Newton's method, starting from `x0`.
    fn solve(self, y: &rug::Float, x0: f64) -> rug::Float {
        let mut x = rug::Float::with_val(EXT_PREC, x0);
        for _ in 0..100 {
            let dx: rug::Float = (self.eval(&x) - y) / self.derivative_prec(&x, EXT_PREC);
            x -= &dx;
            if dx.is_zero() || dx.abs() < rug::Float::with_val(EXT_PREC, x.abs_ref()) >> 128 {
                break;
            }
        }
        x
    }

    /// The `k`-th positive zero of the function (counting from 1) is close to
    /// `(k + zero_offset) * π`.
    fn zero_offset(self) -> f64 {
        match self {
            Self::J0 => -0.25,
            Self::J1 => 0.25,
            Self::Y0 => -0.75,
            Self::Y1 => -0.25,
        }
    }

    /// Returns the value closest to the zero of the function that is the
    /// closest one to `x0`, which must be positive.
    fn zero_near(self, x0: f64) -> f64 {
        // The precision must be enough to tell the distance from the zero to
        // the closest values
        let prec = EXT_PREC + super::exponent(x0).max(0) as u32;
        let pi = rug::Float::with_val(prec, rug::float::Constant::Pi);
        // McMahon's expansion, with b = (k + zero_offset) * π and
        // μ = 4 * order^2: b - (μ - 1) / (8 * b) - ...
        let k = (rug::Float::with_val(prec, x0) / &pi - self.zero_offset())
            .round()
            .max(&rug::Float::with_val(prec, 1));
        let b = (k + self.zero_offset()) * pi;
        let mu = f64::from(4 * self.order().pow(2));
        let mut x: rug::Float = b.clone() - (mu - 1.0) / (8u32 * b);
        // Newton's method
        for _ in 0..100 {
            let dx: rug::Float = self.eval_prec(&x, prec) / self.derivative_prec(&x, prec);
            x -= &dx;
            if dx.is_zero() || dx.abs() < rug::Float::with_val(prec, x.abs_ref()) >> (prec - 32) {
                break;
            }
        }
        x.to_f64()
    }

    /// Returns the `k`-th positive zero of the function (counting from 1),
    /// rounded.
    fn zero(self, k: u32) -> f64 {
        self.zero_near((f64::from(k) + self.zero_offset()) * std::f64::consts::PI)
    }

    /// Returns the extremum of the function that is the closest one to `x0`
    /// (which must be closer to it than to a zero), rounded.
    fn extremum_near(self, x0: f64) -> f64 {
        // Newton's method on the derivative
        let mut x = rug::Float::with_val(EXT_PREC, x0);
        for _ in 0..100 {
            let dx: rug::Float =
                self.derivative_prec(&x, EXT_PREC) / self.derivative2_prec(&x, EXT_PREC);
            x -= &dx;
            if dx.is_zero() || dx.abs() < rug::Float::with_val(EXT_PREC, x.abs_ref()) >> 128 {
                break;
            }
        }
        x.to_f64()
    }

    /// Arguments below which the function is calculated with a polynomial of
    /// `x^2`.
    fn small_max(self) -> f64 {
        match self {
            Self::J0 => 1.0,
            Self::J1 => 2.375,
            Self::Y0 => 0.5,
            Self::Y1 => 1.0,
        }
    }

    /// Arguments from which the function is calculated from its modulus and
    /// phase.
    fn large_min(self) -> f64 {
        match self {
            Self::J0 | Self::Y1 => 47.875,
            Self::J1 | Self::Y0 => 49.375,
        }
    }

    /// Number of zeros below `large_min`.
    fn table_zeros(self) -> u32 {
        match self {
            Self::Y0 => 16,
            _ => 15,
        }
    }
}

fn test_with(func: Func, mut f: impl FnMut(f64)) {
    // Negative arguments are tested by the callers.

    // Test special values
    for x in values::specials() {
        f(x);
    }

    // Test a set of subnormal values, whose results are subnormal with `j1`
    // and overflow with `y1`
    for x in values::subnormals(10_000) {
        f(x);
    }

    // Test a set of normal values at each binade. The values closest to each
    // power of two are included, which covers the thresholds where the
    // calculation changes with tiny and huge arguments (such as 2^-40, 2^35,
    // 2^140 and 2^500).
    for x in values::binades(-1022..=1023, 100) {
        f(x);
    }

    // Test more values with small magnitudes, whose results are calculated
    // with the same polynomials
    for x in values::binades(-45..=-2, 2000) {
        f(x);
    }

    // `rug` is much slower with the functions of the second kind, so they
    // are tested with fewer values
    let scale = match func {
        Func::J0 | Func::J1 => 3,
        Func::Y0 | Func::Y1 => 1,
    };

    // Test values evenly spread over the range of the tables, and beyond
    for x in values::interval(0.0, 52.0, 400_000 * scale) {
        f(x);
    }

    // Test more values where the polynomials of the modulus and the phase
    // are not constant
    for x in values::interval(47.0, 2000.0, 100_000 * scale) {
        f(x);
    }
    for x in values::binades(10..=40, 1500 * scale) {
        f(x);
    }

    // Test the arguments around the limits where the evaluation method
    // changes
    for x in [func.small_max(), func.large_min()] {
        for x in values::around(x, 10_000) {
            f(x);
        }
    }

    // Test the arguments around the limits of the intervals of the tables,
    // which are the midpoints between the approximate positions of the
    // zeros
    for k in 1..=func.table_zeros() {
        let limit = (f64::from(k) + func.zero_offset() + 0.5) * std::f64::consts::PI;
        for x in values::around(limit, 100) {
            f(x);
        }
    }

    // Test the arguments closest to the zeros, whose results are the
    // smallest ones
    for x in near_zeros(func) {
        f(x);
    }

    // Test the arguments closest to the multiples of π/4 at each binade,
    // which are the hardest to reduce: with huge arguments, the zeros and
    // the extrema of the functions are very close to them
    let quarter_pi = |prec| rug::Float::with_val(prec, rug::float::Constant::Pi) / 4;
    for x in values::near_multiples(quarter_pi, 6..=1023) {
        f(x);
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, which are the hardest to round
    for x in near_midpoints(func) {
        f(x);
    }

    // Test integers
    for arg in 1..=20_000 {
        f(f64::from(arg));
    }
}

/// Returns `(1 + m * 2^-52) * 2^e` for some evenly spread `m`: `n` values in
/// each of the binades with exponents in `e`.
fn sparse_binades(e: impl IntoIterator<Item = i16>, n: u32) -> impl Iterator<Item = f64> {
    e.into_iter().flat_map(move |e| {
        values::spread(n, i64::from(e) as u64)
            .map(move |m| fpmath::scalbn(((1u64 << 52) | m) as f64, i32::from(e) - 52))
    })
}

/// Returns the arguments closest to zeros of the function: many arguments
/// around the zeros of the tables, a few around each of the next zeros (up
/// to 2^14), and around some zeros of each of the next binades.
fn near_zeros(func: Func) -> Vec<f64> {
    let mut xs = Vec::new();
    for k in 1..=5000 {
        let n = if k <= func.table_zeros() + 2 { 2000 } else { 2 };
        xs.extend(values::around(func.zero(k), n));
    }
    let more = sparse_binades(14..=70, 40).chain(sparse_binades(71..=1023, 5));
    for x0 in more {
        xs.extend(values::around(func.zero_near(x0), 2));
    }
    xs
}

/// Returns arguments whose results are very close to a midpoint between two
/// consecutive values.
///
/// The functions are flat around their extrema (their results change much
/// less than an ULP between consecutive arguments), so the arguments closest
/// to the inverse of a midpoint have such results (see
/// `values::midpoint_inverse`). One of the extrema of `j0` is zero.
///
/// With small arguments, `2 * j1(x) = x + g(x)` is not flat, but `g` is (see
/// `values::offset_midpoint`).
///
/// Otherwise, such results are searched among the arguments close to some
/// values (see `values::midpoint_scan`), up to where consecutive arguments
/// are too far apart for it.
fn near_midpoints(func: Func) -> Vec<f64> {
    let eval = |x: &rug::Float| func.eval(x);

    let mut xs = Vec::new();

    // The first extrema, which are between the zeros, with distances from
    // 2^-27 to 2^-10 (relative to the extremum)
    let extrema = (1..=40).map(|k| {
        let between = (f64::from(k) + func.zero_offset() + 0.5) * std::f64::consts::PI;
        func.extremum_near(between)
    });
    for xe in extrema {
        for d in sparse_binades(-27..=-10, 3) {
            for x0 in [xe * (1.0 - d), xe * (1.0 + d)] {
                let x = values::midpoint_inverse(x0, eval, |y| func.solve(y, x0));
                xs.extend(values::around(x, 1));
            }
        }
    }

    // The smallest binade from which the function is not flat with small
    // arguments
    let scan_min = match func {
        Func::J0 => {
            for x0 in values::binades(-27..=-4, 300) {
                let x = values::midpoint_inverse(x0, eval, |y| func.solve(y, x0));
                xs.extend(values::around(x, 1));
            }
            -3
        }
        Func::J1 => {
            // g(x) = 2 * j1(x) - x, g'(x) = 2 * j1'(x) - 1
            let twice = |x: &rug::Float| func.eval(x) * 2u32;
            let dg = |x: &rug::Float| func.derivative_prec(x, EXT_PREC) * 2u32 - 1u32;
            for x0 in values::binades(-26..=-2, 300) {
                if let Some(x) = values::offset_midpoint(x0, twice, dg) {
                    xs.extend(values::around(x, 1));
                }
            }
            -1
        }
        // Their magnitudes grow as the argument decreases
        Func::Y0 | Func::Y1 => -1022,
    };

    let small = sparse_binades(scan_min..=5, 4).chain(values::binades(-8..=5, 400));
    let table = values::interval(0.25, 52.0, 5000);
    let large = values::binades(6..=30, 150);
    for x0 in small.chain(table).chain(large) {
        xs.push(values::midpoint_scan(x0, eval, 1 << 12));
    }

    // Close arguments can give the same midpoint
    xs.sort_by(f64::total_cmp);
    xs.dedup();
    xs
}

/// A kind of Bessel function, of any order.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Kind {
    First,
    Second,
}

/// Orders of `jn` and `yn` to test (the callers also test their opposites),
/// where the last ones use the asymptotic expansions of large orders.
const ORDERS: [i32; 18] = [
    2, 3, 4, 5, 6, 7, 8, 10, 13, 20, 33, 64, 100, 250, 601, 1500, 4096, 4099,
];

/// Calculates the Bessel function of the kind `kind` and order `n >= 2` of
/// `x >= 0`, with a relative error of about `2^-prec` (unless the result is
/// extremely close to zero).
///
/// The functions of `rug` take a very long time with large orders and
/// arguments, so this uses the recurrence
///
/// `C[k+1](x) = (2 * k / x) * C[k](x) - C[k-1](x)`
///
/// from the functions of orders 0 and 1, where it is stable: with the
/// functions of the second kind, and with the ones of the first kind when
/// `x >= n`. Its terms nearly cancel when the result is close to zero, so it
/// is calculated with 128 more bits.
///
/// The functions of the first kind use `rug` with `x < n`, where it is fast.
fn eval_order(kind: Kind, n: i32, x: &rug::Float, prec: u32) -> rug::Float {
    assert!(n >= 2);
    if x.is_nan() {
        return rug::Float::with_val(prec, rug::float::Special::Nan);
    }
    if x.is_infinite() || (kind == Kind::First && x.is_zero()) {
        return rug::Float::new(prec);
    }
    if kind == Kind::First && *x < n {
        return rug::Float::with_val(prec, x).jn(n);
    }

    let wp = prec + 128 + 64;
    let xw = rug::Float::with_val(wp, x);
    let (mut a, mut b) = match kind {
        Kind::First => (xw.clone().j0(), xw.clone().j1()),
        Kind::Second => (xw.clone().y0(), xw.clone().y1()),
    };
    for k in 1..n {
        if b.is_infinite() {
            break;
        }
        let c = rug::Float::with_val(wp, &b * (2 * k)) / &xw - &a;
        a = b;
        b = c;
    }
    rug::Float::with_val(prec, b)
}

/// Checks the result of `jn(n, x)` or `yn(n, x)`, for `n >= 2` and `x >= 0`.
///
/// With `x <= n`, the error limit is the one of the other functions, and
/// `max_error` is only updated with those errors. With `x > n`, the functions
/// oscillate, and the error can be larger by 2^-64 times their modulus.
fn check_order_result(kind: Kind, n: i32, x: f64, actual: f64, max_error: &mut f64) {
    let xr = rug::Float::with_val(53, x);
    let eval = |kind| eval_order(kind, n, &xr, RUG_PREC);
    check_order_result_with(kind, n as u32, x, actual, eval, max_error);
}

/// Like `check_order_result`, with the results of `eval` for each kind (only
/// the one of `kind` with `x <= n`).
fn check_order_result_with(
    kind: Kind,
    n: u32,
    x: f64,
    actual: f64,
    eval: impl Fn(Kind) -> rug::Float,
    max_error: &mut f64,
) {
    if x.is_nan() || x <= f64::from(n) || x.is_infinite() {
        check_result((n, x), actual, eval(kind), max_error);
        return;
    }

    // The modulus is sqrt(jn(n, x)^2 + yn(n, x)^2)
    let first = eval(Kind::First);
    let second = eval(Kind::Second);
    let modulus = first.clone().hypot(&second);
    let expected = match kind {
        Kind::First => first,
        Kind::Second => second,
    };

    // In ULPs of the expected result
    let ulp_exp = super::exponent(expected.to_f64()).max(-1022) - 52;
    let extra = (modulus >> (64 + ulp_exp)).to_f64();
    let expected_is_neg = expected.is_sign_negative();
    let err = calc_error_ulp(actual, expected);
    assert!(
        err < ERROR_LIMIT + extra,
        "input = {:?}, error = {err} ULP (limit = {})",
        (n, x),
        ERROR_LIMIT + extra
    );
    // The sign of the smallest results could be wrong
    if extra < 0.25 {
        assert_eq!(
            actual.is_sign_negative(),
            expected_is_neg,
            "input = {:?}, sign mismatch",
            (n, x)
        );
    }
}

/// Returns the arguments closest to the first `count` zeros greater than `n`
/// of the Bessel function of the kind `kind` and order `n >= 2`.
fn order_zeros(kind: Kind, n: i32, count: usize) -> Vec<f64> {
    let is_neg =
        |x: f64| eval_order(kind, n, &rug::Float::with_val(53, x), RUG_PREC).is_sign_negative();

    // The zeros are more than 2 apart, so a step of 1/2 finds all the sign
    // changes
    let mut zeros = Vec::new();
    let mut a = f64::from(n);
    while zeros.len() < count {
        let b = a + 0.5;
        if is_neg(a) != is_neg(b) {
            // Bisection, down to consecutive values
            let (mut lo, mut hi) = (a, b);
            while super::ordinal(hi) - super::ordinal(lo) > 1 {
                let mid = 0.5 * (lo + hi);
                if is_neg(mid) == is_neg(lo) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            zeros.push(lo);
        }
        a = b;
    }
    zeros
}

fn test_order_with(kind: Kind, n: i32, mut f: impl FnMut(f64)) {
    // Negative arguments and orders are tested by the callers.

    let nf = f64::from(n);
    // The recurrences take a time proportional to the order, so the larger
    // ones are tested with fewer values
    let count = |c: u32| (c / (n as u32).max(8) * 8).max(30);

    // Test special values
    for x in values::specials() {
        f(x);
    }

    // Test a set of subnormal values, whose results underflow to zero with
    // `jn` and overflow with `yn`
    for x in values::subnormals(count(300)) {
        f(x);
    }

    // Test a set of normal values at each binade. With tiny values, the
    // results underflow with `jn` and overflow with `yn`. With huge values,
    // they are calculated from their modulus and phase.
    for x in sparse_binades(-1022..=1023, if n <= 64 { 2 } else { 1 }) {
        f(x);
    }

    // Test values evenly spread up to the order, where the function does not
    // oscillate, and well beyond it
    for x in values::interval(0.0, nf, count(8000)) {
        f(x);
    }
    for x in values::interval(nf, 4.0 * nf + 60.0, count(8000)) {
        f(x);
    }
    for x in values::interval(4.0 * nf + 60.0, 40.0 * nf + 2000.0, count(4000)) {
        f(x);
    }

    // Test the arguments around the limits where the evaluation method
    // changes: the order, the argument from which the Debye expansions are
    // used beyond the orders below 1024 (`sqrt(x^2 - n^2) = 646`), and the
    // one from which the functions are calculated from their modulus and
    // phase (n^2 * 2^16)
    let far = fpmath::sqrt(nf * nf + 646.0 * 646.0);
    for x in [nf, far, fpmath::scalbn(nf * nf, 16)] {
        for x in values::around(x, count(200)) {
            f(x);
        }
    }

    // With arguments below the order, the results of `jn` end up being
    // subnormal and then zero as the argument decreases, and the ones of
    // `yn` end up overflowing. Test the arguments around the ones where that
    // starts, and values spread over the results that are close to it.
    let tiny = f64::from_bits(1);
    let eval = |x: f64| eval_order(kind, n, &rug::Float::with_val(53, x), RUG_PREC);
    match kind {
        Kind::First => {
            // The results below this one are rounded to zero
            let half_tiny = rug::Float::with_val(53, tiny) >> 1u32;
            let x_zero = bisect(tiny, nf, |x| eval(x) < half_tiny);
            let x_normal = bisect(tiny, nf, |x| eval(x) < f64::MIN_POSITIVE);
            for x0 in [x_zero, x_normal] {
                for x in values::around(x0, count(200)) {
                    f(x);
                }
            }
            for x in values::interval(x_zero * (1.0 - 1.0 / nf), x_normal, count(400)) {
                f(x);
            }
        }
        Kind::Second => {
            let x_inf = bisect(tiny, nf, |x| eval(x) < -f64::MAX);
            for x in values::around(x_inf, count(200)) {
                f(x);
            }
            // The magnitude of the results decreases by a factor of up to
            // `e` in this interval
            for x in values::interval(x_inf, x_inf * (1.0 + 1.0 / nf), count(400)) {
                f(x);
            }
        }
    }

    // Test the arguments closest to the first zeros, whose results are the
    // smallest ones
    for z in order_zeros(kind, n, 8) {
        for x in values::around(z, 3) {
            f(x);
        }
    }

    // Test results that are very close to a midpoint between two consecutive
    // values, with arguments below the order, where the error limit is the
    // one of the other functions
    let eval = |x: &rug::Float| eval_order(kind, n, x, EXT_PREC);
    for x0 in values::interval(0.25 * nf, nf, count(600)) {
        let x = values::midpoint_scan(x0, eval, 1 << 12);
        if x < nf {
            f(x);
        }
    }
}

/// Returns the first value for which `pred` is false, given that it is true
/// for `lo`, false for `hi`, and it only changes once between them.
fn bisect(mut lo: f64, mut hi: f64, pred: impl Fn(f64) -> bool) -> f64 {
    assert!(pred(lo) && !pred(hi));
    while super::ordinal(hi) - super::ordinal(lo) > 1 {
        let lo_ord = super::ordinal(lo);
        let mid = super::from_ordinal(lo_ord + (super::ordinal(hi) - lo_ord) / 2);
        if pred(mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    hi
}
