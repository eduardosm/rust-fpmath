use crate::traits::{Float, Int as _};

pub(crate) trait Erf: Float {
    fn erf_finite(x: Self) -> Self;

    fn erfc_finite(x: Self) -> Self;

    fn erfcx_finite(x: Self) -> Self;
}

pub(crate) fn erf<F: Erf>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // erf(±0) = ±0
        x
    } else if e == F::MAX_RAW_EXP {
        if x.raw_mant() == F::Raw::ZERO {
            // erf(±inf) = ±1
            F::ONE.copysign(x)
        } else {
            // propagate NaN
            x
        }
    } else {
        // also handles subnormals, erf(x) ~= x * 2 / sqrt(π)
        F::erf_finite(x)
    }
}

pub(crate) fn erfc<F: Erf>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::RawExp::ZERO {
        // x is zero or subnormal
        // erfc(x) ~= 1
        F::ONE
    } else if e == F::MAX_RAW_EXP {
        if x.raw_mant() == F::Raw::ZERO {
            if x.is_sign_negative() {
                // erfc(-inf) = 2
                F::TWO
            } else {
                // erfc(+inf) = 0
                F::ZERO
            }
        } else {
            // propagate NaN
            x
        }
    } else {
        F::erfc_finite(x)
    }
}

pub(crate) fn erfcx<F: Erf>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::RawExp::ZERO {
        // x is zero or subnormal
        // erfcx(x) ~= 1
        F::ONE
    } else if e == F::MAX_RAW_EXP {
        if x.raw_mant() == F::Raw::ZERO {
            if x.is_sign_negative() {
                // erfcx(-inf) = inf
                F::INFINITY
            } else {
                // erfcx(+inf) = 0
                F::ZERO
            }
        } else {
            // propagate NaN
            x
        }
    } else {
        F::erfcx_finite(x)
    }
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::Float;

    fn test_erf<F: Float + FloatMath>() {
        use crate::erf;

        let f = F::parse;

        assert_is_nan!(erf(F::NAN));
        assert_total_eq!(erf(F::INFINITY), F::ONE);
        assert_total_eq!(erf(F::NEG_INFINITY), -F::ONE);
        assert_total_eq!(erf(F::ZERO), F::ZERO);
        assert_total_eq!(erf(-F::ZERO), -F::ZERO);
        assert_total_eq!(erf(f("7")), F::ONE);
        assert_total_eq!(erf(f("-7")), -F::ONE);
        assert_total_eq!(erf(F::LARGEST), F::ONE);
        assert_total_eq!(erf(-F::LARGEST), -F::ONE);
    }

    fn test_erfc<F: Float + FloatMath>() {
        use crate::erfc;

        let f = F::parse;

        assert_is_nan!(erfc(F::NAN));
        assert_total_eq!(erfc(F::INFINITY), F::ZERO);
        assert_total_eq!(erfc(F::NEG_INFINITY), F::TWO);
        assert_total_eq!(erfc(F::ZERO), F::ONE);
        assert_total_eq!(erfc(-F::ZERO), F::ONE);
        assert_total_eq!(erfc(f("28")), F::ZERO);
        assert_total_eq!(erfc(f("-7")), F::TWO);
        assert_total_eq!(erfc(F::LARGEST), F::ZERO);
        assert_total_eq!(erfc(-F::LARGEST), F::TWO);
    }

    fn test_erfcx<F: Float + FloatMath>() {
        use crate::erfcx;

        let f = F::parse;

        assert_is_nan!(erfcx(F::NAN));
        assert_total_eq!(erfcx(F::INFINITY), F::ZERO);
        assert_total_eq!(erfcx(F::NEG_INFINITY), F::INFINITY);
        assert_total_eq!(erfcx(F::ZERO), F::ONE);
        assert_total_eq!(erfcx(-F::ZERO), F::ONE);
        assert_total_eq!(erfcx(f("-27")), F::INFINITY);
        assert_total_eq!(erfcx(-F::LARGEST), F::INFINITY);
        // erfcx(x) ~= 1 / (sqrt(π) * x), which does not underflow to zero
        let tiny = erfcx(F::LARGEST);
        assert!(tiny > F::ZERO && tiny < F::ONE / F::LARGEST);
    }

    #[test]
    fn test_f32() {
        test_erf::<f32>();
        test_erfc::<f32>();
        test_erfcx::<f32>();
    }

    #[test]
    fn test_f64() {
        test_erf::<f64>();
        test_erfc::<f64>();
        test_erfcx::<f64>();
    }
}
