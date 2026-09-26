use crate::traits::{Float, Int as _};

pub(crate) trait InvHyperbolic: Float {
    fn asinh_finite(x: Self) -> Self;

    fn acosh_finite(x: Self) -> Self;

    fn atanh_finite(x: Self) -> Self;
}

pub(crate) fn asinh<F: InvHyperbolic>(x: F) -> F {
    let e = x.raw_exp();
    let m = x.raw_mant();
    if e == F::MAX_RAW_EXP && m != F::Raw::ZERO {
        // asinh(NaN) = NaN
        F::NAN
    } else if (e == F::MAX_RAW_EXP || e == F::RawExp::ZERO) && m == F::Raw::ZERO {
        // asinh(±inf) = ±inf
        // asinh(±0) = ±0
        // asinh(x) = ~x for very small x
        x
    } else {
        F::asinh_finite(x)
    }
}

pub(crate) fn acosh<F: InvHyperbolic>(x: F) -> F {
    let e = x.raw_exp();
    if x < F::ONE {
        // x < 1, acosh(x) is NaN
        F::NAN
    } else if x == F::ONE {
        // asinh(1) = 0
        F::ZERO
    } else if e == F::MAX_RAW_EXP {
        // x is infinity or NaN
        // acosh(x) = x
        x
    } else {
        F::acosh_finite(x)
    }
}

pub(crate) fn atanh<F: InvHyperbolic>(x: F) -> F {
    let absx = x.abs();
    if absx == F::ONE {
        // atanh(±1) = ±inf
        F::INFINITY.copysign(x)
    } else if absx > F::ONE {
        // |x| > 1, return NaN
        F::NAN
    } else {
        F::atanh_finite(x)
    }
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::Float;

    fn test_asinh<F: Float + FloatMath>() {
        use crate::asinh;

        assert_is_nan!(asinh(F::NAN));
        assert_total_eq!(asinh(F::INFINITY), F::INFINITY);
        assert_total_eq!(asinh(F::NEG_INFINITY), F::NEG_INFINITY);
        assert_total_eq!(asinh(F::ZERO), F::ZERO);
        assert_total_eq!(asinh(-F::ZERO), -F::ZERO);
    }

    fn test_acosh<F: Float + FloatMath>() {
        use crate::acosh;

        assert_is_nan!(acosh(F::NAN));
        assert_is_nan!(acosh(F::NEG_INFINITY));
        assert_is_nan!(acosh(-F::ONE));
        assert_is_nan!(acosh(F::ZERO));
        assert_is_nan!(acosh(F::HALF));
        assert_total_eq!(acosh(F::INFINITY), F::INFINITY);
        assert_total_eq!(acosh(F::ONE), F::ZERO);
    }

    fn test_atanh<F: Float + FloatMath>() {
        use crate::atanh;

        let f = F::parse;

        assert_is_nan!(atanh(F::NAN));
        assert_is_nan!(atanh(f("1.5")));
        assert_is_nan!(atanh(f("-1.5")));
        assert_is_nan!(atanh(F::INFINITY));
        assert_is_nan!(atanh(F::NEG_INFINITY));
        assert_total_eq!(atanh(F::ZERO), F::ZERO);
        assert_total_eq!(atanh(-F::ZERO), -F::ZERO);
        assert_total_eq!(atanh(F::ONE), F::INFINITY);
        assert_total_eq!(atanh(-F::ONE), F::NEG_INFINITY);
    }

    #[test]
    fn test_f32() {
        test_asinh::<f32>();
        test_acosh::<f32>();
        test_atanh::<f32>();
    }

    #[test]
    fn test_f64() {
        test_asinh::<f64>();
        test_acosh::<f64>();
        test_atanh::<f64>();
    }
}
