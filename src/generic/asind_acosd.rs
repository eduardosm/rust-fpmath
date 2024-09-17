use super::asin_acos::{acos_inner, asin_inner};
use super::{AsinAcos, RadToDeg};
use crate::traits::Int as _;

pub(crate) fn asind<F: AsinAcos + RadToDeg>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::EXP_OFFSET && x.raw_mant() == F::Raw::ZERO {
        // asind(±1) = ±90
        F::cast_from(90u32).copysign(x)
    } else if e >= F::EXP_OFFSET {
        // NaN or |x| > 1 (including infinity)
        F::NAN
    } else if e <= F::RawExp::from(F::MANT_BITS) {
        // subnormal or zero, asind(x) ~= x * (180/π)
        // also handles asind(-0) = -0
        x * F::rad_to_deg()
    } else {
        let (y_hi, y_lo) = asin_inner(x);
        let (y_hi, y_lo) = F::norm_hi_lo_splitted(y_hi, y_lo);

        y_hi * F::rad_to_deg_hi() + (y_hi * F::rad_to_deg_lo() + y_lo * F::rad_to_deg())
    }
}

pub(crate) fn acosd<F: AsinAcos + RadToDeg>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::EXP_OFFSET && x.raw_mant() == F::Raw::ZERO {
        if x.sign() {
            // acosd(-1) = 180
            F::cast_from(180u32)
        } else {
            // acosd(1) = 0
            F::ZERO
        }
    } else if e >= F::EXP_OFFSET {
        // NaN or |x| > 1 (including infinity)
        F::NAN
    } else if e == F::RawExp::ZERO {
        // subnormal or zero
        // acosd(x) ~= 90
        F::cast_from(90u32)
    } else {
        let (y_hi, y_lo) = acos_inner(x);
        let (y_hi, y_lo) = F::norm_hi_lo_splitted(y_hi, y_lo);

        y_hi * F::rad_to_deg_hi() + (y_hi * F::rad_to_deg_lo() + y_lo * F::rad_to_deg())
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test_asind<F: Float + FloatMath>() {
        use crate::asind;

        let f = F::parse;

        assert_is_nan!(asind(F::NAN));
        assert_is_nan!(asind(f("1.5")));
        assert_is_nan!(asind(f("-1.5")));
        assert_is_nan!(asind(F::INFINITY));
        assert_is_nan!(asind(F::neg_infinity()));
        assert_total_eq!(asind(F::ZERO), F::ZERO);
        assert_total_eq!(asind(-F::ZERO), -F::ZERO);
        assert_total_eq!(asind(F::one()), f("90"));
        assert_total_eq!(asind(-F::one()), f("-90"));
    }

    fn test_acosd<F: Float + FloatMath>() {
        use crate::acosd;

        let f = F::parse;

        assert_is_nan!(acosd(F::NAN));
        assert_is_nan!(acosd(f("1.5")));
        assert_is_nan!(acosd(f("-1.5")));
        assert_is_nan!(acosd(F::INFINITY));
        assert_is_nan!(acosd(F::neg_infinity()));
        assert_total_eq!(acosd(F::ZERO), f("90"));
        assert_total_eq!(acosd(-F::ZERO), f("90"));
        assert_total_eq!(acosd(F::one()), F::ZERO);
        assert_total_eq!(acosd(-F::one()), f("180"));
    }

    #[test]
    fn test_f32() {
        test_asind::<f32>();
        test_acosd::<f32>();
    }

    #[test]
    fn test_soft_f32() {
        test_asind::<crate::SoftF32>();
        test_acosd::<crate::SoftF32>();
    }

    #[test]
    fn test_f64() {
        test_asind::<f64>();
        test_acosd::<f64>();
    }

    #[test]
    fn test_soft_f64() {
        test_asind::<crate::SoftF64>();
        test_acosd::<crate::SoftF64>();
    }
}
