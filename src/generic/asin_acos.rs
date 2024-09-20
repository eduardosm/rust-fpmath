use super::sqrt::two_hi_lo_sqrt_inner;
use crate::traits::{FloatConsts, Int as _, Like};

pub(crate) trait AsinAcos<L = Like<Self>>: FloatConsts {
    fn frac_pi_2_hi() -> Self;
    fn frac_pi_2_lo() -> Self;

    /// Calculates `(asin(x) - x) / x^3`
    fn asin_poly(x2: Self) -> Self;
}

pub(crate) fn asin<F: AsinAcos>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::EXP_OFFSET && x.raw_mant() == F::Raw::ZERO {
        // asin(±1) = ±π/2
        F::frac_pi_2().copysign(x)
    } else if e >= F::EXP_OFFSET {
        // NaN or |x| > 1 (including infinity)
        F::NAN
    } else if e == F::RawExp::ZERO {
        // subnormal or zero, asin(x) ~= x
        // also handles asin(-0) = -0
        x
    } else {
        let (y_hi, y_lo) = asin_inner(x);
        y_hi + y_lo
    }
}

pub(super) fn asin_inner<F: AsinAcos>(x: F) -> (F, F) {
    if x.exponent() < -F::Exp::ONE {
        // |x| < 0.5
        let x2 = x * x;
        let x3 = x2 * x;

        // t1 = asin(x) - x
        let t1 = x3 * F::asin_poly(x2);

        // t2 = asin(x)
        let t2_hi = (t1 + x).purify();
        let t2_lo = (x - t2_hi) + t1;

        (t2_hi, t2_lo)
    } else {
        // |x| >= 0.5
        // |asin(x)| = π/2 - 2 * asin(sqrt((1 - |x|) / 2))

        // y = sqrt((1 - |x|) / 2)
        let y2 = (F::one() - x.abs()) * F::half();
        let (twoy_hi, twoy_lo) = two_hi_lo_sqrt_inner(y2);
        let twoy3 = y2 * twoy_hi;

        // t2 = 2 * (asin(y) - y)
        let t2 = twoy3 * F::asin_poly(y2);

        // t3 = |asin(x)| = π/2 - 2 * asin(y)
        let t3_hi = (F::frac_pi_2_hi() - (twoy_hi + t2)).purify();
        let t3_lo = (F::frac_pi_2_lo() - twoy_lo) + (((F::frac_pi_2_hi() - t3_hi) - twoy_hi) - t2);

        let sgn = F::one().copysign(x);
        (t3_hi * sgn, t3_lo * sgn)
    }
}

pub(super) fn acos_inner<F: AsinAcos>(x: F) -> (F, F) {
    // acos(x) = π/2 - asin(x)
    if x.exponent() < -F::Exp::ONE {
        // |x| < 0.5

        let x2 = x * x;
        let x3 = x2 * x;

        // t1 = asin(x) - x
        let t1 = x3 * F::asin_poly(x2);

        // t2 = asin(x)
        let t2_hi = (F::frac_pi_2_hi() - (t1 + x)).purify();
        let t2_lo = F::frac_pi_2_lo() + (((F::frac_pi_2_hi() - t2_hi) - t1) - x);

        (t2_hi, t2_lo)
    } else {
        // |x| >= 0.5
        // |asin(x)| = π/2 - 2 * asin(sqrt((1 - |x|) / 2))

        // y = sqrt((1 - |x|) / 2)
        let y2 = (F::one() - x.abs()) * F::half();
        let (twoy_hi, twoy_lo) = two_hi_lo_sqrt_inner(y2);
        let twoy3 = y2 * twoy_hi;

        // t1 = 2 * (asin(y) - y)
        let t1 = twoy3 * F::asin_poly(y2);

        // t2 = 2 * asin(y)
        let t2_hi = (t1 + twoy_hi).purify();
        let t2_lo = twoy_lo + ((twoy_hi - t2_hi) + t1);

        if x > F::ZERO {
            // acos(x) = π/2 - |asin(x)|
            //         = π/2 - (π/2 - 2 * asin(y))
            //         = 2 * asin(y)
            (t2_hi, t2_lo)
        } else {
            // acos(x) = π/2 + |asin(x)|
            //         = π/2 + (π/2 - 2 * asin(y))
            //         = π - 2 * asin(y)
            let pi_hi = F::frac_pi_2_hi() * F::two();
            let pi_lo = F::frac_pi_2_lo() * F::two();

            let t3_hi = (pi_hi - t2_hi).purify();
            let t3_lo = (pi_lo - t2_lo) + ((pi_hi - t3_hi) - t2_hi);

            (t3_hi, t3_lo)
        }
    }
}

pub(crate) fn acos<F: AsinAcos>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::EXP_OFFSET && x.raw_mant() == F::Raw::ZERO {
        if x.sign() {
            // acos(-1) = π
            F::pi()
        } else {
            // acos(1) = 0
            F::ZERO
        }
    } else if e >= F::EXP_OFFSET {
        // NaN or |x| > 1 (including infinity)
        F::NAN
    } else if e == F::RawExp::ZERO {
        // subnormal or zero
        // acos(x) ~= π/2
        F::frac_pi_2()
    } else {
        let (y_hi, y_lo) = acos_inner(x);
        y_hi + y_lo
    }
}

#[cfg(test)]
mod tests {
    use super::AsinAcos;
    use crate::FloatMath;

    fn test_asin<F: AsinAcos + FloatMath>() {
        use crate::asin;

        let f = F::parse;

        assert_is_nan!(asin(F::NAN));
        assert_is_nan!(asin(f("1.5")));
        assert_is_nan!(asin(f("-1.5")));
        assert_is_nan!(asin(F::INFINITY));
        assert_is_nan!(asin(F::neg_infinity()));
        assert_total_eq!(asin(F::ZERO), F::ZERO);
        assert_total_eq!(asin(-F::ZERO), -F::ZERO);
        assert_total_eq!(asin(F::one()), F::frac_pi_2());
        assert_total_eq!(asin(-F::one()), -F::frac_pi_2());
    }

    fn test_acos<F: AsinAcos + FloatMath>() {
        use crate::acos;

        let f = F::parse;

        assert_is_nan!(acos(F::NAN));
        assert_is_nan!(acos(f("1.5")));
        assert_is_nan!(acos(f("-1.5")));
        assert_is_nan!(acos(F::INFINITY));
        assert_is_nan!(acos(F::neg_infinity()));
        assert_total_eq!(acos(F::ZERO), F::frac_pi_2());
        assert_total_eq!(acos(-F::ZERO), F::frac_pi_2());
        assert_total_eq!(acos(F::one()), F::ZERO);
        assert_total_eq!(acos(-F::one()), F::pi());
    }

    #[test]
    fn test_f32() {
        test_asin::<f32>();
        test_acos::<f32>();
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f32() {
        test_asin::<crate::SoftF32>();
        test_acos::<crate::SoftF32>();
    }

    #[test]
    fn test_f64() {
        test_asin::<f64>();
        test_acos::<f64>();
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f64() {
        test_asin::<crate::SoftF64>();
        test_acos::<crate::SoftF64>();
    }
}
