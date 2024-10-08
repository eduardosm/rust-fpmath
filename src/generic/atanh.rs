use super::log::log_hi_lo_inner;
use super::Log;
use crate::double::SemiDouble;
use crate::traits::Int as _;

pub(crate) fn atanh<F: Log>(x: F) -> F {
    let e = x.raw_exp();
    let absx = x.abs();
    if (e == F::MAX_RAW_EXP && x.raw_mant() != F::Raw::ZERO)
        || e <= (F::EXP_OFFSET - F::RawExp::from(F::MANT_BITS))
    {
        // propagate NaN
        // or
        // very small, includes subnormal and zero
        // atanh(x) ~= x
        // also handles atanh(-0) = -0
        x
    } else if absx == F::one() {
        // atanh(±1) = ±inf
        F::INFINITY.copysign(x)
    } else if x.abs() > F::one() {
        // |x| > 1, return NaN
        F::NAN
    } else {
        atanh_inner(x)
    }
}

fn atanh_inner<F: Log>(x: F) -> F {
    // t1 = 2 * x / (1 - x)
    let t1 = SemiDouble::new(F::two() * x) / SemiDouble::new_qsub11(F::one(), x);

    // t2 = (1 + x) / (1 - x) = t1 + 1
    let t2 = t1 + F::one();
    let t2 = t2.to_norm();

    // atanh(x) = 0.5 * log((1 + x) / (1 - x))
    F::half() * log_hi_lo_inner(t2.hi(), t2.lo())
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test<F: Float + FloatMath>() {
        use crate::atanh;

        let f = F::parse;

        assert_is_nan!(atanh(F::NAN));
        assert_is_nan!(atanh(f("1.5")));
        assert_is_nan!(atanh(f("-1.5")));
        assert_is_nan!(atanh(F::INFINITY));
        assert_is_nan!(atanh(F::neg_infinity()));
        assert_total_eq!(atanh(F::ZERO), F::ZERO);
        assert_total_eq!(atanh(-F::ZERO), -F::ZERO);
        assert_total_eq!(atanh(F::one()), F::INFINITY);
        assert_total_eq!(atanh(-F::one()), F::neg_infinity());
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
