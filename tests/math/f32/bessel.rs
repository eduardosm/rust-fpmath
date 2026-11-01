use super::{MIN_MAX_ERROR, check_result, values};

#[test]
fn test_j0() {
    let mut max_error: f32 = 0.0;
    test_with(Func::J0, |x| {
        let expected = fpmath::j0(f64::from(x));
        let actual = fpmath::j0(x);
        assert_result_eq!(fpmath::j0(-x), actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_j1() {
    let mut max_error: f32 = 0.0;
    test_with(Func::J1, |x| {
        let expected = fpmath::j1(f64::from(x));
        let actual = fpmath::j1(x);
        assert_result_eq!(fpmath::j1(-x), -actual);

        check_result(x, actual, expected, &mut max_error);
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_y0() {
    let mut max_error: f32 = 0.0;
    test_with(Func::Y0, |x| {
        let expected = fpmath::y0(f64::from(x));
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
    let mut max_error: f32 = 0.0;
    test_with(Func::Y1, |x| {
        let expected = fpmath::y1(f64::from(x));
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
    let mut max_error: f32 = 0.0;
    let mut check = |n: i32, x: f32| {
        let expected = fpmath::jn(n, f64::from(x));
        let actual = fpmath::jn(n, x);
        // jn(-n, x) = (-1)^n * jn(n, x) and jn(n, -x) = (-1)^n * jn(n, x)
        let odd = if n % 2 == 0 { actual } else { -actual };
        assert_result_eq!(fpmath::jn(-n, x), odd);
        assert_result_eq!(fpmath::jn(n, -x), odd);
        assert_result_eq!(fpmath::jn(-n, -x), actual);

        check_result((n, x), actual, expected, &mut max_error);
    };
    for n in ORDERS {
        test_order_with(n, |x| check(n, x));
    }
    test_below_order_with(LAST_PLAIN_ORDER, 1000, |x| check(LAST_PLAIN_ORDER, x));
    test_order_with(FIRST_F64_ORDER, |x| check(FIRST_F64_ORDER, x));
    for n in LARGE_ORDERS {
        test_order_with(n, |x| check(n, x));
    }
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

#[test]
fn test_yn() {
    let mut max_error: f32 = 0.0;
    let mut check = |n: i32, x: f32| {
        let expected = fpmath::yn(n, f64::from(x));
        let actual = fpmath::yn(n, x);
        // yn(-n, x) = (-1)^n * yn(n, x)
        let odd = if n % 2 == 0 { actual } else { -actual };
        assert_result_eq!(fpmath::yn(-n, x), odd);
        if x > 0.0 {
            assert!(fpmath::yn(n, -x).is_nan());
        }

        check_result((n, x), actual, expected, &mut max_error);
    };
    for n in ORDERS {
        test_order_with(n, |x| check(n, x));
    }
    test_below_order_with(LAST_PLAIN_ORDER, 1000, |x| check(LAST_PLAIN_ORDER, x));
    test_order_with(FIRST_F64_ORDER, |x| check(FIRST_F64_ORDER, x));
    for n in LARGE_ORDERS {
        test_order_with(n, |x| check(n, x));
    }
    eprintln!("max error = {max_error}");
    assert!(max_error > MIN_MAX_ERROR);
}

/// A Bessel function of the first (`J`) or second (`Y`) kind, of order 0 or 1.
#[derive(Copy, Clone, PartialEq, Eq)]
enum Func {
    J0,
    J1,
    Y0,
    Y1,
}

impl Func {
    fn order(self) -> u32 {
        match self {
            Self::J0 | Self::Y0 => 0,
            Self::J1 | Self::Y1 => 1,
        }
    }

    /// Calculates the function in `f64`.
    fn eval(self, x: f64) -> f64 {
        match self {
            Self::J0 => fpmath::j0(x),
            Self::J1 => fpmath::j1(x),
            Self::Y0 => fpmath::y0(x),
            Self::Y1 => fpmath::y1(x),
        }
    }

    /// Calculates the derivative of the function in `f64`, with
    /// `C0'(x) = -C1(x)` and `C1'(x) = C0(x) - C1(x) / x`.
    fn derivative(self, x: f64) -> f64 {
        match self {
            Self::J0 => -fpmath::j1(x),
            Self::J1 => fpmath::j0(x) - fpmath::j1(x) / x,
            Self::Y0 => -fpmath::y1(x),
            Self::Y1 => fpmath::y0(x) - fpmath::y1(x) / x,
        }
    }

    /// Solves `func(x) = y` in `f64` with Newton's method, starting from
    /// `x0`.
    fn solve(self, y: f64, x0: f64) -> f64 {
        let mut x = x0;
        for _ in 0..100 {
            let dx = (self.eval(x) - y) / self.derivative(x);
            x -= dx;
            if dx.abs() <= x.abs() * f64::EPSILON {
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
    fn zero_near(self, x0: f64) -> f32 {
        // McMahon's expansion, with b = (k + zero_offset) * π and
        // μ = 4 * order^2: b - (μ - 1) / (8 * b) - ...
        let pi = std::f64::consts::PI;
        let k = fpmath::round(x0 / pi - self.zero_offset()).max(1.0);
        let b = (k + self.zero_offset()) * pi;
        let mu = f64::from(4 * self.order().pow(2));
        self.solve(0.0, b - (mu - 1.0) / (8.0 * b)) as f32
    }

    /// Arguments below which the function is calculated with a polynomial of
    /// `x^2`.
    fn small_max(self) -> f32 {
        match self {
            Self::J0 => 1.0,
            Self::J1 => 2.375,
            Self::Y0 => 0.5,
            Self::Y1 => 1.0,
        }
    }

    /// Arguments from which the function is calculated from its modulus and
    /// phase.
    fn large_min(self) -> f32 {
        match self {
            Self::J0 | Self::Y1 => 47.875,
            Self::J1 | Self::Y0 => 49.375,
        }
    }
}

fn test_with(func: Func, mut f: impl FnMut(f32)) {
    // Negative arguments are tested by the callers.

    // Test special values
    for x in values::specials() {
        f(x);
    }

    // Test a set of subnormal values, whose results are subnormal with `j1`
    // and overflow with `y1`
    for x in values::subnormals(100_000) {
        f(x);
    }

    // Test a set of normal values at each binade
    for x in values::binades(-126..=127, 1000) {
        f(x);
    }

    // Test more values with small magnitudes, whose results are calculated
    // with the same polynomials (only with their first terms below 2^-27,
    // the values closest to each power of two are included)
    for x in values::binades(-30..=-1, 100_000) {
        f(x);
    }

    // Exhaustive test of all the arguments from 1 up to where the functions
    // are calculated from their modulus and phase, and beyond. With `y0`,
    // from the end of its polynomial of small arguments, which includes its
    // first zero (about 0.89, the others are larger than 1 for all the
    // functions).
    for x in values::interval_full(func.small_max().min(1.0), 64.0) {
        f(x);
    }

    // Test more values with larger magnitudes
    for x in values::interval(64.0, 1.0e6, 5_000_000) {
        f(x);
    }
    for x in values::binades(6..=127, 20_000) {
        f(x);
    }

    // Test the arguments around the limits where the evaluation method
    // changes
    for x in [func.small_max(), func.large_min()] {
        for x in values::around(x, 100_000) {
            f(x);
        }
    }

    // Test the arguments closest to the zeros beyond the exhaustive range,
    // whose results are the smallest ones: all of them up to 2^16, and some
    // of them up to 2^24 (where the spacing of the values is larger than the
    // one of the zeros)
    for k in 20..=21_000 {
        let z = func.zero_near((f64::from(k) + func.zero_offset()) * std::f64::consts::PI);
        for x in values::around(z, 2) {
            f(x);
        }
    }
    for x0 in values::binades(16..=23, 5000) {
        for x in values::around(func.zero_near(f64::from(x0)), 1) {
            f(x);
        }
    }

    // Test the arguments closest to the multiples of π/4 at each binade,
    // which are the hardest to reduce: with huge arguments, the zeros and
    // the extrema of the functions are very close to them
    let quarter_pi = |prec| rug::Float::with_val(prec, rug::float::Constant::Pi) / 4;
    for x in values::near_multiples(quarter_pi, 6..=127) {
        f(x);
    }

    // Test results that are very close to a midpoint between two consecutive
    // values with small arguments, which are the hardest to round
    for x in near_midpoints(func) {
        f(x);
    }

    // Test integers
    for arg in 1..=200_000u32 {
        f(arg as f32);
    }
}

/// Returns small arguments whose results are very close to a midpoint between
/// two consecutive values.
///
/// `j0` is flat around zero (its result changes much less than an ULP between
/// consecutive arguments), so the arguments closest to the inverse of a
/// midpoint have such results (see `values::midpoint_inverse`).
///
/// `2 * j1(x) = x + g(x)` is not flat, but `g` is (see
/// `values::offset_midpoint`).
fn near_midpoints(func: Func) -> Vec<f32> {
    let eval = |x: f64| func.eval(x);

    let mut xs = Vec::new();
    match func {
        Func::J0 => {
            for x0 in values::binades(-13..=-3, 2000) {
                let x = values::midpoint_inverse(x0, eval, |y| func.solve(y, f64::from(x0)));
                xs.extend(values::around(x, 1));
            }
        }
        Func::J1 => {
            // g(x) = 2 * j1(x) - x, g'(x) = 2 * j1'(x) - 1
            let twice = |x: f64| 2.0 * func.eval(x);
            let dg = |x: f64| 2.0 * func.derivative(x) - 1.0;
            for x0 in values::binades(-12..=-2, 2000) {
                if let Some(x) = values::offset_midpoint(x0, twice, dg) {
                    xs.extend(values::around(x, 1));
                }
            }
        }
        Func::Y0 | Func::Y1 => {}
    }

    // Close arguments can give the same midpoint
    xs.sort_by(f32::total_cmp);
    xs.dedup();
    xs
}

/// Orders of `jn` and `yn` to test (the callers also test their opposites).
const ORDERS: [i32; 9] = [2, 3, 4, 5, 8, 21, 64, 300, 700];

/// The largest order whose functions have their own recurrences for
/// arguments below the order. They are calculated with plain `f64`
/// arithmetic, whose error grows with the order, so this one has the largest
/// errors. It is only tested with those arguments, because of the time that
/// the `f64` functions take with it.
const LAST_PLAIN_ORDER: i32 = 1023;

/// The first order whose functions use the `f64` ones for all the arguments.
/// From it, those use asymptotic expansions, which take a time that does not
/// depend on the order, or know beforehand that the result is zero or
/// infinite as `f32`.
const FIRST_F64_ORDER: i32 = LAST_PLAIN_ORDER + 1;

/// Orders above `FIRST_F64_ORDER`, up to the largest one.
const LARGE_ORDERS: [i32; 3] = [16_384, 1_000_003, i32::MAX];

fn test_order_with(n: i32, mut f: impl FnMut(f32)) {
    // Negative arguments and orders are tested by the callers.

    // Test special values
    for x in values::specials() {
        f(x);
    }

    // Test a set of subnormal values, whose results underflow to zero with
    // `jn` and overflow with `yn`
    for x in values::subnormals(order_count(n, 10_000)) {
        f(x);
    }

    // Test a set of normal values at each binade. With small values, the
    // results underflow with `jn` (to subnormal values or zero) and overflow
    // with `yn`.
    for x in values::binades(-126..=127, order_count(n, 1000).min(200)) {
        f(x);
    }

    test_above_order_with(n, &mut f);
    test_below_order_with(n, 1000, f);
}

/// Tests arguments from the order `n`, where the functions oscillate.
fn test_above_order_with(n: i32, mut f: impl FnMut(f32)) {
    let nf = n as f32;

    // Test values evenly spread well beyond the order
    for x in values::interval(nf, 40.0 * nf + 2000.0, order_count(n, 400_000)) {
        f(x);
    }

    // Test the arguments around the one from which the functions are
    // calculated from their modulus and phase (n^2 * 2^16)
    for x in values::around(fpmath::scalbn(nf * nf, 16), order_count(n, 1000)) {
        f(x);
    }
}

/// Tests arguments up to the order `n`, where the functions do not oscillate.
///
/// `min_count` is the least number of values to test between the arguments
/// from which the results are not zero or infinite and the order.
fn test_below_order_with(n: i32, min_count: u32, mut f: impl FnMut(f32)) {
    let nf = n as f32;
    let tiny = f32::from_bits(1);

    // Test values evenly spread up to the order
    for x in values::interval(0.0, nf, order_count(n, 400_000)) {
        f(x);
    }

    // Test the arguments around the order, where the evaluation method
    // changes
    for x in values::around(nf, order_count(n, 1000)) {
        f(x);
    }

    // With smaller arguments, `jn` ends up underflowing to zero and `yn`
    // overflowing, which is most of the arguments below the order with a
    // large one. Test the arguments around the ones where that starts, and
    // values evenly spread from there up to the order.
    let jn_is_zero = |x: f32| fpmath::jn(n, f64::from(x)) < fpmath::scalbn(1.0, -150);
    let yn_is_inf = |x: f32| fpmath::yn(n, f64::from(x)) < -f64::from(f32::MAX);
    for x0 in [bisect(tiny, nf, jn_is_zero), bisect(tiny, nf, yn_is_inf)] {
        for x in values::around(x0, 100) {
            f(x);
        }
        for x in values::interval(x0, nf, order_count(n, 100_000).max(min_count)) {
            f(x);
        }
    }
}

/// Returns how many values to test out of `c` with the order `n`: up to
/// `LAST_PLAIN_ORDER`, the functions take a time proportional to the order
/// (and so do the `f64` ones that the tests compare with), so the larger ones
/// are tested with fewer values. From `FIRST_F64_ORDER`, the time does not
/// depend on the order (also with the `f64` functions), so they are tested
/// with all of them.
fn order_count(n: i32, c: u32) -> u32 {
    if n >= FIRST_F64_ORDER {
        return c;
    }
    (c / (n as u32).max(8) * 8).max(100)
}

/// Returns the first value for which `pred` is false, given that it is true
/// for `lo`, false for `hi`, and it only changes once between them.
fn bisect(mut lo: f32, mut hi: f32, pred: impl Fn(f32) -> bool) -> f32 {
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
