use crate::traits::{Float, Int as _};

pub(crate) trait Log: Float {
    fn ln_finite(x: Self, edelta: Self::Exp) -> Self;

    fn ln_1p_finite(x: Self) -> Self;

    fn log2_finite(x: Self, edelta: Self::Exp) -> Self;

    fn log10_finite(x: Self, edelta: Self::Exp) -> Self;
}

pub(crate) fn ln<F: Log>(x: F) -> F {
    let (y, edelta) = x.normalize_arg();
    let yexp = y.raw_exp();
    if yexp == F::RawExp::ZERO {
        // ln(±0) = -inf
        F::NEG_INFINITY
    } else if y.is_sign_negative() {
        // x < 0, ln(x) = NaN
        F::NAN
    } else if yexp == F::MAX_RAW_EXP {
        // propagate infinity or NaN
        y
    } else {
        F::ln_finite(y, edelta)
    }
}

pub(crate) fn ln_1p<F: Log>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::RawExp::ZERO {
        // subnormal or zero, log(1 + x) ~= x
        // also handles log(1 + (-0)) = -0
        x
    } else if x == -F::ONE {
        // x = -1, log(1 + x) = -inf
        F::NEG_INFINITY
    } else if x < -F::ONE {
        // x < -1, log(1 + x) = NaN
        F::NAN
    } else if e == F::MAX_RAW_EXP {
        // propagate infinity or NaN
        x
    } else {
        F::ln_1p_finite(x)
    }
}

pub(crate) fn log2<F: Log>(x: F) -> F {
    let (y, edelta) = x.normalize_arg();
    let yexp = y.raw_exp();
    if yexp == F::RawExp::ZERO {
        // log2(±0) = -inf
        F::NEG_INFINITY
    } else if y.is_sign_negative() {
        // x < 0, log2(x) = NaN
        F::NAN
    } else if yexp == F::MAX_RAW_EXP {
        if y.raw_mant() == F::Raw::ZERO {
            // log2(inf) = inf
            F::INFINITY
        } else {
            // NaN, propagate
            y
        }
    } else {
        F::log2_finite(y, edelta)
    }
}

pub(crate) fn log10<F: Log>(x: F) -> F {
    let (y, edelta) = x.normalize_arg();
    let yexp = y.raw_exp();
    if yexp == F::RawExp::ZERO {
        // log10(±0) = -inf
        F::NEG_INFINITY
    } else if y.is_sign_negative() {
        // x < 0, log10(x) = NaN
        F::NAN
    } else if yexp == F::MAX_RAW_EXP {
        if y.raw_mant() == F::Raw::ZERO {
            // log10(inf) = inf
            F::INFINITY
        } else {
            // NaN, propagate
            y
        }
    } else {
        F::log10_finite(y, edelta)
    }
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::Float;

    fn test_ln<F: Float + FloatMath>() {
        use crate::ln;

        assert_is_nan!(ln(F::NAN));
        assert_is_nan!(ln(-F::ONE));
        assert_is_nan!(ln(F::NEG_INFINITY));
        assert_total_eq!(ln(F::ZERO), F::NEG_INFINITY);
        assert_total_eq!(ln(-F::ZERO), F::NEG_INFINITY);
        assert_total_eq!(ln(F::INFINITY), F::INFINITY);
    }

    fn test_ln_1p<F: Float + FloatMath>() {
        use crate::ln_1p;

        assert_is_nan!(ln_1p(F::NAN));
        assert_is_nan!(ln_1p(-(F::ONE + F::HALF)));
        assert_is_nan!(ln_1p(F::NEG_INFINITY));
        assert_total_eq!(ln_1p(-F::ONE), F::NEG_INFINITY);
        assert_total_eq!(ln_1p(-F::ZERO), -F::ZERO);
        assert_total_eq!(ln_1p(F::ZERO), F::ZERO);
        assert_total_eq!(ln_1p(F::INFINITY), F::INFINITY);
    }

    fn test_log2<F: Float + FloatMath>() {
        use crate::log2;

        assert_is_nan!(log2(F::NAN));
        assert_is_nan!(log2(-F::ONE));
        assert_is_nan!(log2(F::NEG_INFINITY));
        assert_total_eq!(log2(F::ZERO), F::NEG_INFINITY);
        assert_total_eq!(log2(-F::ZERO), F::NEG_INFINITY);
        assert_total_eq!(log2(F::INFINITY), F::INFINITY);
    }

    fn test_log10<F: Float + FloatMath>() {
        use crate::log10;

        assert_is_nan!(log10(F::NAN));
        assert_is_nan!(log10(-F::ONE));
        assert_is_nan!(log10(F::NEG_INFINITY));
        assert_total_eq!(log10(F::ZERO), F::NEG_INFINITY);
        assert_total_eq!(log10(-F::ZERO), F::NEG_INFINITY);
        assert_total_eq!(log10(F::INFINITY), F::INFINITY);
    }

    #[test]
    fn test_f32() {
        test_ln::<f32>();
        test_ln_1p::<f32>();
        test_log2::<f32>();
        test_log10::<f32>();
    }

    #[test]
    fn test_f64() {
        test_ln::<f64>();
        test_ln_1p::<f64>();
        test_log2::<f64>();
        test_log10::<f64>();
    }
}
