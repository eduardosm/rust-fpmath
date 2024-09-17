use crate::traits::{CastInto as _, Float, Int as _, Like};

pub(crate) trait Log<L = Like<Self>>: Float {
    fn sqrt_2() -> Self;
    fn ln_2_hi() -> Self;
    fn ln_2_lo() -> Self;
    fn frac_2_3_hi() -> Self;
    fn frac_2_3_lo() -> Self;
    fn frac_4_10_hi() -> Self;
    fn frac_4_10_lo() -> Self;

    /// Calculates `(log(1 + x) - log(1 - x) - 2 * x) / x`
    ///
    /// `-0.1716 < x < 0.1716`
    fn log_special_poly(x: Self) -> Self;

    /// Calculates `(log(1 + x) - log(1 - x) - 2 * x - (2/3) * x^3 - 0.4 * x^5)
    /// / x`
    ///
    /// `-0.1716 < x < 0.1716`
    fn log_special_poly_ex(x2: Self) -> Self;
}

/// Returns `(k, r)` as needed by `log_inner`
///
/// * `x * 2^edelta = 2^k * (1 + r)`
/// * `sqrt(2) / 2 <= r + 1 <= sqrt(2)`
pub(super) fn log_split<F: Log>(x: F, edelta: F::Exp) -> (F::Exp, F) {
    // Split x * 2^edelta = 2^k * m
    //  - k is an integer
    //  - 1 <= m < 2
    let k = x.exponent() + edelta;
    // Take the mantissa from x and set the expontent to 0
    let m = x.set_exp(F::Exp::ZERO);

    // reduce 1 <= m < 2 into sqrt(2) / 2 <= 1 + r <= sqrt(2)
    if m > F::sqrt_2() {
        (k + F::Exp::ONE, m.set_exp(-F::Exp::ONE) - F::one())
    } else {
        (k, m - F::one())
    }
}

/// Calculates `log(x * 2^edelta)`
///
/// `x` must be normal and positive.
pub(super) fn log_inner<F: Log>(x: F, edelta: F::Exp) -> F {
    // Algorithm based on one used by the msun math library:
    //  * log(1 + r) = p * s + 2 * s
    //  * s = r / (2 + r)
    //  * p = (log(1 + s) - log(1 - s) - 2 * s) / s

    // Split x * 2^edelta = 2^k * (1 + r)
    //  - k is an integer
    //  - sqrt(2) / 2 <= 1 + r < sqrt(2)
    let (k, r) = log_split(x, edelta);

    // s = r / (2 + r)
    // So, log(1 + r) = log(1 + s) - log(1 - s)
    let s = r / (F::two() + r);

    // p = (log(1 + s) - log(1 - s) - 2 * s) / s
    let p = F::log_special_poly(s);

    // t1 = k * log(2)
    let kf: F = k.cast_into();
    let t1_hi = kf * F::ln_2_hi();
    let t1_lo = kf * F::ln_2_lo();

    // log(x) = log(1 + r) + k * log(2) = log(1 + r) + t1
    // where log(1 + r) = p * s + 2 * s
    //                  = r - s * (r - p)
    //                  = r - (0.5 * r^2 - s * (0.5 * r^2 + p))
    let hr2 = F::half() * r * r;
    (((s * (hr2 + p) + t1_lo) - hr2) + r) + t1_hi
}

pub(crate) fn log<F: Log>(x: F) -> F {
    let (y, edelta) = x.normalize_arg();
    let yexp = y.raw_exp();
    if yexp == F::RawExp::ZERO {
        // log(±0) = -inf
        F::neg_infinity()
    } else if y.sign() {
        // x < 0, log(x) = NaN
        F::NAN
    } else if yexp == F::MAX_RAW_EXP {
        // propagate infinity or NaN
        y
    } else {
        log_inner(y, edelta)
    }
}

// Calculates log(x_hi + x_lo)
pub(super) fn log_hi_lo_inner<F: Log>(x_hi: F, x_lo: F) -> F {
    // Algorithm based on one used by the msun math library:
    //  * log(1 + r) = p * s + 2 * s
    //  * s = r / (2 + r)
    //  * p = (log(1 + s) - log(1 - s) - 2 * s) / s

    // Split x_hi = 2^k * (1 + r)
    //  - k is an integer
    //  - sqrt(2) / 2 <= 1 + r < sqrt(2)
    //  - e is an error term
    let (k, r) = log_split(x_hi, F::Exp::ZERO);

    // Calculate a correction term to handle x_lo:
    // log(x_hi + x_lo) = log(x_hi) + c
    // c = log(x_hi + x_lo) - log(x_hi) =
    //   = log((x_hi + x_lo) / x_hi) =
    //   = log(1 + x_lo / x_hi) ~= x_lo / x_hi
    let c = x_lo / x_hi;

    // s = r / (2 + r)
    // So, log(1 + r) = log(1 + s) - log(1 - s)
    let s = r / (F::two() + r);

    // p = (log(1 + s) - log(1 - s) - 2 * s) / s
    let p = F::log_special_poly(s);

    // t1 = k * log(2)
    let kf: F = k.cast_into();
    let t1_hi = kf * F::ln_2_hi();
    let t1_lo = kf * F::ln_2_lo();

    // log(x) = log(x_hi) + c
    //        = log(1 + r) + k * log(2) + c
    //        = log(1 + r) + t1 + c
    // log(1 + r) = p * s + 2 * s
    //            = r - s * (r - p)
    //            = r - (0.5 * r^2 - s * (0.5 * r^2 + p))
    let hr2 = F::half() * r * r;
    (((s * (hr2 + p) + (t1_lo + c)) - hr2) + r) + t1_hi
}

fn log_1p_inner<F: Log>(x: F) -> F {
    // Calculate xp1 + e = 1 + x, where e is an
    // error term to handle rounding in 1 + x.
    let xp1 = (F::one() + x).purify();
    let e = if x > F::one() {
        (x - xp1) + F::one()
    } else {
        (F::one() - xp1) + x
    };

    // Calculate log(1 + x) = log(xp1 + e)
    log_hi_lo_inner(xp1, e)
}

pub(crate) fn log_1p<F: Log>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::RawExp::ZERO {
        // subnormal or zero, log(1 + x) ~= x
        // also handles log(1 + (-0)) = -0
        x
    } else if x == -F::one() {
        // x = -1, log(1 + x) = -inf
        F::neg_infinity()
    } else if x < -F::one() {
        // x < -1, log(1 + x) = NaN
        F::NAN
    } else if e == F::MAX_RAW_EXP {
        // propagate infinity or NaN
        x
    } else {
        log_1p_inner(x)
    }
}

/// Calculates `(log(1 + x) - log(1 - x) - 2 * x) / x`
///
/// `-0.1716 < x < 0.1716`
pub(super) fn hi_lo_log_special_poly<F: Log>(x2_hi: F, x2_lo: F) -> (F, F) {
    // p1 = (p - 2/3 * x^2 - 0.4 * x^4) / x^4
    let x2 = x2_hi + x2_lo;
    let p0 = F::log_special_poly_ex(x2);

    // p1 = (p - 2/3 * x^2) / x^4
    let p1_hi = (p0 + F::frac_4_10_hi()).purify();
    let p1_lo = ((F::frac_4_10_hi() - p1_hi) + p0) + F::frac_4_10_lo();
    let (p1_hi, p1_lo) = F::norm_hi_lo_splitted(p1_hi, p1_lo);

    // p2 = (p - 2/3 * x^2) / x^2
    let p2_hi = p1_hi * x2_hi;
    let p2_lo = p1_hi * x2_lo + p1_lo * x2_hi + p1_lo * x2_lo;

    // p3 = p / x^2
    let p3_hi = (p2_hi + F::frac_2_3_hi()).purify();
    let p3_lo = (((F::frac_2_3_hi() - p3_hi) + p2_hi) + F::frac_2_3_lo()) + p2_lo;
    let (p3_hi, p3_lo) = F::norm_hi_lo_splitted(p3_hi, p3_lo);

    // p = (log(1 + x) - log(1 - x) - 2 * x) / x
    let p_hi = p3_hi * x2_hi;
    let p_lo = p3_hi * x2_lo + p3_lo * x2_hi + p3_lo * x2_lo;

    (p_hi, p_lo)
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test_log<F: Float + FloatMath>() {
        use crate::log;

        assert_is_nan!(log(F::NAN));
        assert_is_nan!(log(-F::one()));
        assert_is_nan!(log(F::neg_infinity()));
        assert_total_eq!(log(F::ZERO), F::neg_infinity());
        assert_total_eq!(log(-F::ZERO), F::neg_infinity());
        assert_total_eq!(log(F::INFINITY), F::INFINITY);
    }

    fn test_log_1p<F: Float + FloatMath>() {
        use crate::log_1p;

        assert_is_nan!(log_1p(F::NAN));
        assert_is_nan!(log_1p(-(F::one() + F::half())));
        assert_is_nan!(log_1p(F::neg_infinity()));
        assert_total_eq!(log_1p(-F::one()), F::neg_infinity());
        assert_total_eq!(log_1p(-F::ZERO), -F::ZERO);
        assert_total_eq!(log_1p(F::ZERO), F::ZERO);
        assert_total_eq!(log_1p(F::INFINITY), F::INFINITY);
    }

    #[test]
    fn test_f32() {
        test_log::<f32>();
        test_log_1p::<f32>();
    }

    #[test]
    fn test_soft_f32() {
        test_log::<crate::SoftF32>();
        test_log_1p::<crate::SoftF32>();
    }

    #[test]
    fn test_f64() {
        test_log::<f64>();
        test_log_1p::<f64>();
    }

    #[test]
    fn test_soft_f64() {
        test_log::<crate::SoftF64>();
        test_log_1p::<crate::SoftF64>();
    }
}
