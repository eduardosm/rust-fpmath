use super::is_pos_normal;
use crate::traits::{Float, Int as _};

pub(crate) trait Log: Float {
    fn ln_finite(x: Self, edelta: Self::Exp) -> Self;

    fn ln_1p_finite(x: Self) -> Self;

    fn log2_finite(x: Self, edelta: Self::Exp) -> Self;

    fn log2_1p_finite(x: Self) -> Self;

    fn log10_finite(x: Self, edelta: Self::Exp) -> Self;

    fn log10_1p_finite(x: Self) -> Self;
}

pub(crate) fn ln<F: Log>(x: F) -> F {
    if is_pos_normal(x) {
        // fast path for the most common case
        F::ln_finite(x, F::Exp::ZERO)
    } else {
        log_special(x, F::ln_finite)
    }
}

pub(crate) fn ln_1p<F: Log>(x: F) -> F {
    let e = x.raw_exp();
    if x > -F::ONE && e != F::RawExp::ZERO && e != F::MAX_RAW_EXP {
        // fast path for the most common case (also excludes NaN)
        F::ln_1p_finite(x)
    } else if e == F::RawExp::ZERO {
        // subnormal or zero, log(1 + x) ~= x
        // also handles log(1 + (-0)) = -0
        x
    } else if x == -F::ONE {
        // x = -1, log(1 + x) = -inf
        F::NEG_INFINITY
    } else if x < -F::ONE {
        // x < -1, log(1 + x) = NaN
        F::NAN
    } else {
        // propagate infinity or NaN
        x
    }
}

pub(crate) fn log2<F: Log>(x: F) -> F {
    if is_pos_normal(x) {
        // fast path for the most common case
        F::log2_finite(x, F::Exp::ZERO)
    } else {
        log_special(x, F::log2_finite)
    }
}

pub(crate) fn log2_1p<F: Log>(x: F) -> F {
    let e = x.raw_exp();
    if x > -F::ONE && e != F::RawExp::ZERO && e != F::MAX_RAW_EXP {
        // fast path for the most common case (also excludes NaN)
        F::log2_1p_finite(x)
    } else if e == F::RawExp::ZERO {
        if x.raw_mant() == F::Raw::ZERO {
            // log2(1 + (±0)) = ±0
            x
        } else {
            // subnormal, log2(1 + x) ~= x * log2(e)
            F::log2_1p_finite(x)
        }
    } else if x == -F::ONE {
        // x = -1, log2(1 + x) = -inf
        F::NEG_INFINITY
    } else if x < -F::ONE {
        // x < -1, log2(1 + x) = NaN
        F::NAN
    } else {
        // propagate infinity or NaN
        x
    }
}

pub(crate) fn log10<F: Log>(x: F) -> F {
    if is_pos_normal(x) {
        // fast path for the most common case
        F::log10_finite(x, F::Exp::ZERO)
    } else {
        log_special(x, F::log10_finite)
    }
}

pub(crate) fn log10_1p<F: Log>(x: F) -> F {
    let e = x.raw_exp();
    if x > -F::ONE && e != F::RawExp::ZERO && e != F::MAX_RAW_EXP {
        // fast path for the most common case (also excludes NaN)
        F::log10_1p_finite(x)
    } else if e == F::RawExp::ZERO {
        if x.raw_mant() == F::Raw::ZERO {
            // log10(1 + (±0)) = ±0
            x
        } else {
            // subnormal, log10(1 + x) ~= x * log10(e)
            F::log10_1p_finite(x)
        }
    } else if x == -F::ONE {
        // x = -1, log10(1 + x) = -inf
        F::NEG_INFINITY
    } else if x < -F::ONE {
        // x < -1, log10(1 + x) = NaN
        F::NAN
    } else {
        // propagate infinity or NaN
        x
    }
}

/// Calculates `ln`, `log2` or `log10` (with `finite`) of `x` when it is not
/// positive and normal.
#[inline]
fn log_special<F: Float>(x: F, finite: impl FnOnce(F, F::Exp) -> F) -> F {
    let (y, edelta) = x.normalize_arg();
    let yexp = y.raw_exp();
    if yexp == F::RawExp::ZERO {
        // log(±0) = -inf
        F::NEG_INFINITY
    } else if y.is_sign_negative() {
        // x < 0, log(x) = NaN
        F::NAN
    } else if yexp == F::MAX_RAW_EXP {
        // propagate infinity or NaN
        y
    } else {
        // positive subnormal
        finite(y, edelta)
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

    fn test_log2_1p<F: Float + FloatMath>() {
        use crate::log2_1p;

        let f = F::parse;

        assert_is_nan!(log2_1p(F::NAN));
        assert_is_nan!(log2_1p(-(F::ONE + F::HALF)));
        assert_is_nan!(log2_1p(F::NEG_INFINITY));
        assert_total_eq!(log2_1p(-F::ONE), F::NEG_INFINITY);
        assert_total_eq!(log2_1p(-F::ZERO), -F::ZERO);
        assert_total_eq!(log2_1p(F::ZERO), F::ZERO);
        assert_total_eq!(log2_1p(F::INFINITY), F::INFINITY);
        assert_total_eq!(log2_1p(F::ONE), F::ONE);
        assert_total_eq!(log2_1p(f("3")), F::TWO);
        assert_total_eq!(log2_1p(f("65535")), f("16"));
        assert_total_eq!(log2_1p(-F::HALF), -F::ONE);
        assert_total_eq!(log2_1p(f("-0.75")), -F::TWO);
        assert_total_eq!(log2_1p(f("-0.875")), f("-3"));
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

    fn test_log10_1p<F: Float + FloatMath>() {
        use crate::log10_1p;

        let f = F::parse;

        assert_is_nan!(log10_1p(F::NAN));
        assert_is_nan!(log10_1p(-(F::ONE + F::HALF)));
        assert_is_nan!(log10_1p(F::NEG_INFINITY));
        assert_total_eq!(log10_1p(-F::ONE), F::NEG_INFINITY);
        assert_total_eq!(log10_1p(-F::ZERO), -F::ZERO);
        assert_total_eq!(log10_1p(F::ZERO), F::ZERO);
        assert_total_eq!(log10_1p(F::INFINITY), F::INFINITY);
        assert_total_eq!(log10_1p(f("9")), F::ONE);
        assert_total_eq!(log10_1p(f("99")), F::TWO);
        assert_total_eq!(log10_1p(f("999999")), f("6"));
    }

    #[test]
    fn test_f32() {
        test_ln::<f32>();
        test_ln_1p::<f32>();
        test_log2::<f32>();
        test_log2_1p::<f32>();
        test_log10::<f32>();
        test_log10_1p::<f32>();
    }

    #[test]
    fn test_f64() {
        test_ln::<f64>();
        test_ln_1p::<f64>();
        test_log2::<f64>();
        test_log2_1p::<f64>();
        test_log10::<f64>();
        test_log10_1p::<f64>();
    }
}
