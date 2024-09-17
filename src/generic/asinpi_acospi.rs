use super::asin_acos::{acos_inner, asin_inner};
use super::{AsinAcos, DivPi};
use crate::traits::{CastFrom as _, Int as _};

pub(crate) fn asinpi<F: AsinAcos + DivPi>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::EXP_OFFSET && x.raw_mant() == F::Raw::ZERO {
        // asinpi(±1) = ±0.5
        F::half().copysign(x)
    } else if e >= F::EXP_OFFSET {
        // NaN or |x| > 1 (including infinity)
        F::NAN
    } else if e == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // asinpi(±0) = ±0
        x
    } else if e <= F::RawExp::from(F::MANT_BITS) {
        // very small, asinpi(x) ~= x / π

        // scale temporarily to avoid temporary subnormal numbers
        let logscale = F::Exp::TWO * F::Exp::cast_from(F::MANT_BITS);
        let scale = F::exp2i_fast(logscale);
        let descale = F::exp2i_fast(-logscale);

        let (x_hi, x_lo) = (x * scale).split_hi_lo();
        let y_hi = x_hi * F::frac_1_pi_hi();
        let y_lo = x_hi * F::frac_1_pi_lo() + x_lo * F::frac_1_pi();
        (y_hi + y_lo) * descale
    } else {
        let (y_hi, y_lo) = asin_inner(x);
        let (y_hi, y_lo) = F::norm_hi_lo_splitted(y_hi, y_lo);

        y_hi * F::frac_1_pi_hi() + (y_hi * F::frac_1_pi_lo() + y_lo * F::frac_1_pi())
    }
}

pub(crate) fn acospi<F: AsinAcos + DivPi>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::EXP_OFFSET && x.raw_mant() == F::Raw::ZERO {
        if x.sign() {
            // acospi(-1) = 1
            F::one()
        } else {
            // acospi(1) = 0
            F::ZERO
        }
    } else if e >= F::EXP_OFFSET {
        // NaN or |x| > 1 (including infinity)
        F::NAN
    } else if e == F::RawExp::ZERO {
        // subnormal or zero
        // acospi(x) ~= 0.5
        F::half()
    } else {
        let (y_hi, y_lo) = acos_inner(x);
        let (y_hi, y_lo) = F::norm_hi_lo_splitted(y_hi, y_lo);

        y_hi * F::frac_1_pi_hi() + (y_hi * F::frac_1_pi_lo() + y_lo * F::frac_1_pi())
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test_asinpi<F: Float + FloatMath>() {
        use crate::asinpi;

        let f = F::parse;

        assert_is_nan!(asinpi(F::NAN));
        assert_is_nan!(asinpi(f("1.5")));
        assert_is_nan!(asinpi(f("-1.5")));
        assert_is_nan!(asinpi(F::INFINITY));
        assert_is_nan!(asinpi(F::neg_infinity()));
        assert_total_eq!(asinpi(F::ZERO), F::ZERO);
        assert_total_eq!(asinpi(-F::ZERO), -F::ZERO);
        assert_total_eq!(asinpi(F::one()), F::half());
        assert_total_eq!(asinpi(-F::one()), -F::half());
    }

    fn test_acospi<F: Float + FloatMath>() {
        use crate::acospi;

        let f = F::parse;

        assert_is_nan!(acospi(F::NAN));
        assert_is_nan!(acospi(f("1.5")));
        assert_is_nan!(acospi(f("-1.5")));
        assert_is_nan!(acospi(F::INFINITY));
        assert_is_nan!(acospi(F::neg_infinity()));
        assert_total_eq!(acospi(F::ZERO), F::half());
        assert_total_eq!(acospi(-F::ZERO), F::half());
        assert_total_eq!(acospi(F::one()), F::ZERO);
        assert_total_eq!(acospi(-F::one()), F::one());
    }

    #[test]
    fn test_f32() {
        test_asinpi::<f32>();
        test_acospi::<f32>();
    }

    #[test]
    fn test_soft_f32() {
        test_asinpi::<crate::SoftF32>();
        test_acospi::<crate::SoftF32>();
    }

    #[test]
    fn test_f64() {
        test_asinpi::<f64>();
        test_acospi::<f64>();
    }

    #[test]
    fn test_soft_f64() {
        test_asinpi::<crate::SoftF64>();
        test_acospi::<crate::SoftF64>();
    }
}
