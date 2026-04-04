use crate::double::{DenormDouble, NormDouble, SemiDouble};
use crate::traits::{CastInto as _, Float, Int as _};

pub(crate) trait Ln: Float {
    fn sqrt_2() -> Self;
    fn ln_2_hi() -> Self;
    fn ln_2_lo() -> Self;
    fn frac_2_3_ex() -> NormDouble<Self>;
    fn frac_4_10_ex() -> NormDouble<Self>;

    /// Calculates `(ln(1 + x) - ln(1 - x) - 2 * x) / x`
    ///
    /// `-0.1716 < x < 0.1716`
    fn ln_special_poly(x: Self) -> Self;

    /// Calculates `(ln(1 + x) - ln(1 - x) - 2 * x - (2/3) * x^3 - 0.4 * x^5) / x`
    ///
    /// `-0.1716 < x < 0.1716`
    fn ln_special_poly_ex(x2: Self) -> Self;
}

pub(crate) fn ln<F: Ln>(x: F) -> F {
    let (y, edelta) = x.normalize_arg();
    let yexp = y.raw_exp();
    if yexp == F::RawExp::ZERO {
        // ln(±0) = -inf
        F::neg_infinity()
    } else if y.sign() {
        // x < 0, ln(x) = NaN
        F::NAN
    } else if yexp == F::MAX_RAW_EXP {
        // propagate infinity or NaN
        y
    } else {
        ln_inner(y, edelta)
    }
}

pub(crate) fn ln_1p<F: Ln>(x: F) -> F {
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
        ln_1p_inner(x)
    }
}

/// Calculates `ln(x * 2^edelta)`
///
/// `x` must be normal and positive.
pub(super) fn ln_inner<F: Ln>(x: F, edelta: F::Exp) -> F {
    // Algorithm based on one used by the msun math library:
    //  * ln(1 + r) = p * s + 2 * s
    //  * s = r / (2 + r)
    //  * p = (ln(1 + s) - ln(1 - s) - 2 * s) / s

    // Split x * 2^edelta = 2^k * (1 + r)
    //  - k is an integer
    //  - sqrt(2) / 2 <= 1 + r < sqrt(2)
    let (k, r) = ln_split(x, edelta);

    // s = r / (2 + r)
    // So, ln(1 + r) = ln(1 + s) - ln(1 - s)
    let s = r / (F::two() + r);

    // p = (ln(1 + s) - ln(1 - s) - 2 * s) / s
    let p = F::ln_special_poly(s);

    // t1 = k * ln(2)
    let kf: F = k.cast_into();
    let t1_hi = kf * F::ln_2_hi();
    let t1_lo = kf * F::ln_2_lo();

    // ln(x) = ln(1 + r) + k * ln(2) = ln(1 + r) + t1
    // where ln(1 + r) = p * s + 2 * s
    //                 = r - s * (r - p)
    //                 = r - (0.5 * r^2 - s * (0.5 * r^2 + p))
    let hr2 = F::half() * r * r;
    (((s * (hr2 + p) + t1_lo) - hr2) + r) + t1_hi
}

fn ln_1p_inner<F: Ln>(x: F) -> F {
    // Calculate xp1 + e = 1 + x, where e is an
    // error term to handle rounding in 1 + x.
    let xp1 = (F::one() + x).purify();
    let e = if x > F::one() {
        (x - xp1) + F::one()
    } else {
        (F::one() - xp1) + x
    };

    // Calculate ln(1 + x) = ln(xp1 + e)
    ln_hi_lo_inner(xp1, e)
}

/// Returns `(k, r)` as needed by `ln_inner`
///
/// * `x * 2^edelta = 2^k * (1 + r)`
/// * `sqrt(2) / 2 <= r + 1 <= sqrt(2)`
pub(super) fn ln_split<F: Ln>(x: F, edelta: F::Exp) -> (F::Exp, F) {
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

// Calculates ln(x_hi + x_lo)
pub(super) fn ln_hi_lo_inner<F: Ln>(x_hi: F, x_lo: F) -> F {
    // Algorithm based on one used by the msun math library:
    //  * ln(1 + r) = p * s + 2 * s
    //  * s = r / (2 + r)
    //  * p = (ln(1 + s) - ln(1 - s) - 2 * s) / s

    // Split x_hi = 2^k * (1 + r)
    //  - k is an integer
    //  - sqrt(2) / 2 <= 1 + r < sqrt(2)
    //  - e is an error term
    let (k, r) = ln_split(x_hi, F::Exp::ZERO);

    // Calculate a correction term to handle x_lo:
    // ln(x_hi + x_lo) = ln(x_hi) + c
    // c = ln(x_hi + x_lo) - ln(x_hi) =
    //   = ln((x_hi + x_lo) / x_hi) =
    //   = ln(1 + x_lo / x_hi) ~= x_lo / x_hi
    let c = x_lo / x_hi;

    // s = r / (2 + r)
    // So, ln(1 + r) = ln(1 + s) - ln(1 - s)
    let s = r / (F::two() + r);

    // p = (ln(1 + s) - ln(1 - s) - 2 * s) / s
    let p = F::ln_special_poly(s);

    // t1 = k * ln(2) + c
    let kf: F = k.cast_into();
    let t1 = DenormDouble::new(F::ln_2_hi(), F::ln_2_lo())
        .pmul1(kf)
        .ladd(c);

    // ln(x) = ln(x_hi) + c
    //       = ln(1 + r) + k * ln(2) + c
    //       = ln(1 + r) + t1
    // ln(1 + r) = p * s + 2 * s
    //           = r - s * (r - p)
    //           = r - (0.5 * r^2 - s * (0.5 * r^2 + p))
    let hr2 = F::half() * r * r;
    (((s * (hr2 + p) + t1.lo()) - hr2) + r) + t1.hi()
}

/// Calculates ln(x * 2^edelta)
pub(super) fn hi_lo_ln_inner<F: Ln>(x: F, edelta: F::Exp) -> DenormDouble<F> {
    // Algorithm based on one used by the msun math library:
    //  * ln(1 + r) = p * s + 2 * s
    //  * s = r / (2 + r)
    //  * p = (ln(1 + s) - ln(1 - s) - 2 * s) / s

    // Split x * 2^edelta = 2^k * (1 + r)
    //  - k is an integer
    //  - sqrt(2) / 2 <= 1 + r < sqrt(2)
    let (k, r) = ln_split(x, edelta);

    // rp2 = 2 + r
    let rp2 = SemiDouble::new_qadd11(F::two(), r);

    // s = r / (2 + r)
    let s = (SemiDouble::new(r) / rp2).to_semi();
    let s2 = s.square().to_semi();

    // p = (ln(1 + s) - ln(1 - s) - 2 * s) / s
    let p = hi_lo_ln_special_poly(s2).to_semi();

    // t1 = k * ln(2)
    let kf: F = k.cast_into();
    let t1 = DenormDouble::new(F::ln_2_hi(), F::ln_2_lo()).pmul1(kf);

    // t2 = ln(1 + r) = p * s + 2 * s
    let ps = p * s;
    let twos = s.pmul1(F::two());
    let t2 = twos.to_denorm().qadd2(ps);

    // ln(2^k * (1 + r)) = t1 + t2
    t1.qadd2(t2)
}

/// Calculates ln((x_hi + x_lo) * 2^edelta)
pub(super) fn hi_lo_ln_hi_lo_inner<F: Ln>(x: NormDouble<F>, edelta: F::Exp) -> DenormDouble<F> {
    // Algorithm based on one used by the msun math library:
    //  * ln(1 + r) = p * s + 2 * s
    //  * s = r / (2 + r)
    //  * p = (ln(1 + s) - ln(1 - s) - 2 * s) / s

    // Split x_hi * 2^edelta = 2^k * (1 + r)
    //  - k is an integer
    //  - sqrt(2) / 2 <= 1 + r < sqrt(2)
    let (k, r) = ln_split(x.hi(), edelta);

    // Calculate a correction term to handle x_lo:
    // ln(x_hi + x_lo) = ln(x_hi) + c
    // c = ln(x_hi + x_lo) - ln(x_hi) =
    //   = ln((x_hi + x_lo) / x_hi) =
    //   = ln(1 + x_lo / x_hi) ~= x_lo / x_hi
    let c = x.lo() / x.hi();

    // rp2 = 2 + r
    let rp2 = SemiDouble::new_qadd11(F::two(), r);

    // s = r / (2 + r)
    let s = (SemiDouble::new(r) / rp2).to_semi();
    let s2 = s.square().to_semi();

    // p = (ln(1 + s) - ln(1 - s) - 2 * s) / s
    let p = hi_lo_ln_special_poly(s2).to_semi();

    // t1 = k * ln(2) + c
    let kf: F = k.cast_into();
    let t1 = DenormDouble::new(F::ln_2_hi(), F::ln_2_lo())
        .pmul1(kf)
        .ladd(c);

    // t2 = ln(1 + r) = p * s + 2 * s
    let ps = p * s;
    let twos = s.pmul1(F::two());
    let t2 = twos.to_denorm().qadd2(ps);

    // ln(2^k * (1 + r)) + c = t1 + t2
    t1.qadd2(t2)
}

/// Calculates `(ln(1 + x) - ln(1 - x) - 2 * x) / x`
///
/// `-0.1716 < x < 0.1716`
fn hi_lo_ln_special_poly<F: Ln>(x2: SemiDouble<F>) -> DenormDouble<F> {
    // p0 = (p - 2/3 * x^2 - 0.4 * x^4) / x^4
    let p0 = F::ln_special_poly_ex(x2.to_single());

    // p1 = (p - 2/3 * x^2) / x^4 = p0 + 0.4
    let p1 = F::frac_4_10_ex().to_denorm().qadd1(p0).to_semi();

    // p2 = (p - 2/3 * x^2) / x^2 = p1 * x2
    let p2 = p1 * x2;

    // p3 = p / x^2 = p2 + 2/3
    let p3 = F::frac_2_3_ex().to_denorm().qadd2(p2).to_semi();

    // (log(1 + x) - log(1 - x) - 2 * x) / x = p3 * x2
    p3 * x2
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::Float;

    fn test_ln<F: Float + FloatMath>() {
        use crate::ln;

        assert_is_nan!(ln(F::NAN));
        assert_is_nan!(ln(-F::one()));
        assert_is_nan!(ln(F::neg_infinity()));
        assert_total_eq!(ln(F::ZERO), F::neg_infinity());
        assert_total_eq!(ln(-F::ZERO), F::neg_infinity());
        assert_total_eq!(ln(F::INFINITY), F::INFINITY);
    }

    fn test_ln_1p<F: Float + FloatMath>() {
        use crate::ln_1p;

        assert_is_nan!(ln_1p(F::NAN));
        assert_is_nan!(ln_1p(-(F::one() + F::half())));
        assert_is_nan!(ln_1p(F::neg_infinity()));
        assert_total_eq!(ln_1p(-F::one()), F::neg_infinity());
        assert_total_eq!(ln_1p(-F::ZERO), -F::ZERO);
        assert_total_eq!(ln_1p(F::ZERO), F::ZERO);
        assert_total_eq!(ln_1p(F::INFINITY), F::INFINITY);
    }

    #[test]
    fn test_f32() {
        test_ln::<f32>();
        test_ln_1p::<f32>();
    }

    #[test]
    fn test_f64() {
        test_ln::<f64>();
        test_ln_1p::<f64>();
    }
}
