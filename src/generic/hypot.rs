use super::sqrt::hi_lo_sqrt_hi_lo_inner;
use crate::traits::{Float, Int as _};

pub(crate) fn hypot<F: Float>(x: F, y: F) -> F {
    let xexp = x.raw_exp();
    let yexp = y.raw_exp();
    if xexp == F::MAX_RAW_EXP || yexp == F::MAX_RAW_EXP {
        if xexp == F::MAX_RAW_EXP && x.raw_mant() == F::Raw::ZERO
            || yexp == F::MAX_RAW_EXP && y.raw_mant() == F::Raw::ZERO
        {
            // x or y is inf
            F::INFINITY
        } else {
            // x or y is NaN
            F::NAN
        }
    } else {
        // min = min(|x|, |y|)
        // max = max(|x|, |y|)
        let absx = x.abs();
        let absy = y.abs();
        let (min, max) = if absx < absy {
            (absx, absy)
        } else {
            (absy, absx)
        };

        let maxexp = max.exponent();
        let logscale = maxexp.clamp(F::MIN_NORMAL_EXP, -F::MIN_NORMAL_EXP);
        let scale = F::exp2i_fast(-logscale);
        let descale = F::exp2i_fast(logscale);

        let smin = min * scale;
        let smax = max * scale;

        if smax.raw_exp() == F::RawExp::ZERO {
            F::ZERO
        } else {
            // hypot(x, y) = hypot(min, max)
            //             = hypot(min * scale, max * scale) / scale
            let (smin_hi, smin_lo) = smin.split_hi_lo();
            let (smax_hi, smax_lo) = smax.split_hi_lo();

            let smin2_hi = smin_hi * smin_hi;
            let smin2_lo = F::two() * smin_hi * smin_lo + smin_lo * smin_lo;
            let (smin2_hi, smin2_lo) = F::norm_hi_lo_full(smin2_hi, smin2_lo);

            let smax2_hi = smax_hi * smax_hi;
            let smax2_lo = F::two() * smax_hi * smax_lo + smax_lo * smax_lo;
            let (smax2_hi, smax2_lo) = F::norm_hi_lo_full(smax2_hi, smax2_lo);

            // sum = (min * scale)^2 + (max * scale)^2
            let sum_hi = (smin2_hi + smax2_hi).purify();
            let sum_lo = (((smax2_hi - sum_hi) + smin2_hi) + smax2_lo) + smin2_lo;

            // z = sqrt((min * scale)^2 + (max * scale)^2)
            //   = hypot(min * scale, max * scale)
            let (z_hi, z_lo) = hi_lo_sqrt_hi_lo_inner(sum_hi, sum_lo);

            // hypot(x, y) = hypot(min * scale, max * scale) / scale
            (z_hi + z_lo) * descale
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test<F: Float + FloatMath>() {
        use crate::hypot;

        let f = F::parse;

        assert_is_nan!(hypot(F::NAN, F::NAN));
        assert_is_nan!(hypot(F::NAN, F::ZERO));
        assert_is_nan!(hypot(F::NAN, F::one()));
        assert_is_nan!(hypot(F::ZERO, F::NAN));
        assert_is_nan!(hypot(F::one(), F::NAN));
        assert_total_eq!(hypot(F::INFINITY, F::INFINITY), F::INFINITY);
        assert_total_eq!(hypot(F::INFINITY, F::neg_infinity()), F::INFINITY);
        assert_total_eq!(hypot(F::neg_infinity(), F::INFINITY), F::INFINITY);
        assert_total_eq!(hypot(F::neg_infinity(), F::neg_infinity()), F::INFINITY);
        assert_total_eq!(hypot(F::INFINITY, F::NAN), F::INFINITY);
        assert_total_eq!(hypot(F::neg_infinity(), F::NAN), F::INFINITY);
        assert_total_eq!(hypot(F::NAN, F::INFINITY), F::INFINITY);
        assert_total_eq!(hypot(F::NAN, F::neg_infinity()), F::INFINITY);
        assert_total_eq!(hypot(F::INFINITY, F::ZERO), F::INFINITY);
        assert_total_eq!(hypot(F::INFINITY, F::one()), F::INFINITY);
        assert_total_eq!(hypot(F::neg_infinity(), F::ZERO), F::INFINITY);
        assert_total_eq!(hypot(F::neg_infinity(), F::one()), F::INFINITY);
        assert_total_eq!(hypot(F::ZERO, F::INFINITY), F::INFINITY);
        assert_total_eq!(hypot(F::one(), F::INFINITY), F::INFINITY);
        assert_total_eq!(hypot(F::ZERO, F::neg_infinity()), F::INFINITY);
        assert_total_eq!(hypot(F::one(), F::neg_infinity()), F::INFINITY);
        assert_total_eq!(hypot(F::ZERO, F::ZERO), F::ZERO);
        assert_total_eq!(hypot(-F::ZERO, -F::ZERO), F::ZERO);
        assert_total_eq!(hypot(F::one(), F::ZERO), F::one());
        assert_total_eq!(hypot(f("3"), F::ZERO), f("3"));
        assert_total_eq!(hypot(F::ZERO, F::one()), F::one());
        assert_total_eq!(hypot(F::ZERO, f("3")), f("3"));
        assert_total_eq!(hypot(f("3"), f("4")), f("5"));
    }

    #[test]
    fn test_f32() {
        test::<f32>();
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f32() {
        test::<crate::SoftF32>();
    }

    #[test]
    fn test_f64() {
        test::<f64>();
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f64() {
        test::<crate::SoftF64>();
    }
}
