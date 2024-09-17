use super::{reduce_90_deg, tan::tan_inner, Reduce90Deg, Tan};

pub(crate) fn tand<F: Reduce90Deg + Tan>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        // tand(inf or NaN) = NaN
        F::NAN
    } else if e <= F::RawExp::from(F::MANT_BITS) {
        // very small, includes subnormal and zero
        // tand(x) ~= x * (π/180)
        // also handles tand(-0) = -0
        x * F::deg_to_rad()
    } else {
        let (n, y_hi, y_lo) = reduce_90_deg(x);
        let inv = (n & 1) != 0;
        if inv && y_hi == F::ZERO {
            F::INFINITY.copysign(x)
        } else {
            tan_inner(y_hi, y_lo, inv)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test<F: Float + FloatMath>() {
        use crate::tand;

        assert_is_nan!(tand(F::NAN));
        assert_is_nan!(tand(F::INFINITY));
        assert_is_nan!(tand(F::neg_infinity()));
        assert_total_eq!(tand(F::ZERO), F::ZERO);
        assert_total_eq!(tand(-F::ZERO), -F::ZERO);
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
