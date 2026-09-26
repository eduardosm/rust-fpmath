use crate::traits::{Float, Int as _};

pub(crate) trait Exp: Float {
    fn exp_finite(x: Self) -> Self;

    fn exp_m1_finite(x: Self) -> Self;

    fn exp2_finite(x: Self) -> Self;

    fn exp10_finite(x: Self) -> Self;
}

pub(crate) fn exp<F: Exp>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::RawExp::ZERO {
        // x is zero or subnormal
        // exp(x) ~= 1
        F::ONE
    } else if e == F::MAX_RAW_EXP {
        if x.raw_mant() == F::Raw::ZERO {
            if x.is_sign_negative() {
                // exp(-inf) = 0
                F::ZERO
            } else {
                // exp(+inf) = +inf
                F::INFINITY
            }
        } else {
            // propagate NaN
            x
        }
    } else {
        F::exp_finite(x)
    }
}

pub(crate) fn exp_m1<F: Exp>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::RawExp::ZERO {
        // x is zero or subnormal
        // exp_m1(x) ~= x
        x
    } else if e == F::MAX_RAW_EXP {
        if x.raw_mant() == F::Raw::ZERO {
            if x.is_sign_negative() {
                // exp_m1(-inf) = -1
                -F::ONE
            } else {
                // exp_m1(+inf) = +inf
                F::INFINITY
            }
        } else {
            // propagate NaN
            x
        }
    } else {
        F::exp_m1_finite(x)
    }
}

pub(crate) fn exp2<F: Exp>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::RawExp::ZERO {
        // x is zero or subnormal
        // exp2(x) ~= 1
        F::ONE
    } else if e == F::MAX_RAW_EXP {
        if x.raw_mant() == F::Raw::ZERO {
            if x.is_sign_negative() {
                // exp2(-inf) = 0
                F::ZERO
            } else {
                // exp2(+inf) = +inf
                F::INFINITY
            }
        } else {
            // propagate NaN
            x
        }
    } else {
        F::exp2_finite(x)
    }
}

pub(crate) fn exp10<F: Exp>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::RawExp::ZERO {
        // x is zero or subnormal
        // exp10(x) ~= 1
        F::ONE
    } else if e == F::MAX_RAW_EXP {
        if x.raw_mant() == F::Raw::ZERO {
            if x.is_sign_negative() {
                // exp10(-inf) = 0
                F::ZERO
            } else {
                // exp10(+inf) = +inf
                F::INFINITY
            }
        } else {
            // propagate NaN
            x
        }
    } else {
        F::exp10_finite(x)
    }
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::Float;

    fn test_exp<F: Float + FloatMath>(lo_th: F, hi_th: F) {
        use crate::exp;

        assert_is_nan!(exp(F::NAN));
        assert_total_eq!(exp(F::INFINITY), F::INFINITY);
        assert_total_eq!(exp(F::NEG_INFINITY), F::ZERO);
        assert_total_eq!(exp(F::ZERO), F::ONE);
        assert_total_eq!(exp(-F::ZERO), F::ONE);
        assert_total_eq!(exp(lo_th), F::ZERO);
        assert_total_eq!(exp(lo_th - F::ONE), F::ZERO);
        assert_total_eq!(exp(lo_th - F::TWO), F::ZERO);
        assert_total_eq!(exp(hi_th), F::INFINITY);
        assert_total_eq!(exp(hi_th + F::ONE), F::INFINITY);
        assert_total_eq!(exp(hi_th + F::TWO), F::INFINITY);
    }

    fn test_exp_m1<F: Float + FloatMath>(lo_th: F, hi_th: F) {
        use crate::exp_m1;

        assert_is_nan!(exp_m1(F::NAN));
        assert_total_eq!(exp_m1(F::INFINITY), F::INFINITY);
        assert_total_eq!(exp_m1(F::NEG_INFINITY), -F::ONE);
        assert_total_eq!(exp_m1(F::ZERO), F::ZERO);
        assert_total_eq!(exp_m1(-F::ZERO), -F::ZERO);
        assert_total_eq!(exp_m1(lo_th), -F::ONE);
        assert_total_eq!(exp_m1(lo_th - F::ONE), -F::ONE);
        assert_total_eq!(exp_m1(lo_th - F::TWO), -F::ONE);
        assert_total_eq!(exp_m1(hi_th), F::INFINITY);
        assert_total_eq!(exp_m1(hi_th + F::ONE), F::INFINITY);
        assert_total_eq!(exp_m1(hi_th + F::TWO), F::INFINITY);
    }

    fn test_exp2<F: Float + FloatMath>(lo_th: F, hi_th: F) {
        use crate::exp2;

        let f = F::parse;

        assert_is_nan!(exp2(F::NAN));
        assert_total_eq!(exp2(F::INFINITY), F::INFINITY);
        assert_total_eq!(exp2(F::NEG_INFINITY), F::ZERO);
        assert_total_eq!(exp2(F::ZERO), F::ONE);
        assert_total_eq!(exp2(-F::ZERO), F::ONE);
        assert_total_eq!(exp2(F::ONE), F::TWO);
        assert_total_eq!(exp2(F::TWO), f("4"));
        assert_total_eq!(exp2(f("32")), f("4294967296"));
        assert_total_eq!(exp2(-F::ONE), F::HALF);
        assert_total_eq!(exp2(-F::TWO), f("0.25"));
        assert_total_eq!(exp2(f("-3")), f("0.125"));
        assert_total_eq!(exp2(f("-4")), f("0.0625"));
        assert_total_eq!(exp2(lo_th), F::ZERO);
        assert_total_eq!(exp2(lo_th - F::ONE), F::ZERO);
        assert_total_eq!(exp2(lo_th - F::TWO), F::ZERO);
        assert_total_eq!(exp2(hi_th), F::INFINITY);
        assert_total_eq!(exp2(hi_th + F::ONE), F::INFINITY);
        assert_total_eq!(exp2(hi_th + F::TWO), F::INFINITY);
    }

    fn test_exp10<F: Float + FloatMath>(lo_th: F, hi_th: F) {
        use crate::exp10;

        let f = F::parse;

        assert_is_nan!(exp10(F::NAN));
        assert_total_eq!(exp10(F::INFINITY), F::INFINITY);
        assert_total_eq!(exp10(F::NEG_INFINITY), F::ZERO);
        assert_total_eq!(exp10(F::ZERO), F::ONE);
        assert_total_eq!(exp10(-F::ZERO), F::ONE);
        assert_total_eq!(exp10(F::ONE), f("10"));
        assert_total_eq!(exp10(F::TWO), f("100"));
        assert_total_eq!(exp10(lo_th), F::ZERO);
        assert_total_eq!(exp10(lo_th - F::ONE), F::ZERO);
        assert_total_eq!(exp10(lo_th - F::TWO), F::ZERO);
        assert_total_eq!(exp10(hi_th), F::INFINITY);
        assert_total_eq!(exp10(hi_th + F::ONE), F::INFINITY);
        assert_total_eq!(exp10(hi_th + F::TWO), F::INFINITY);
    }

    #[test]
    fn test_f32() {
        test_exp::<f32>(-103.99, 88.9);
        test_exp_m1::<f32>(-87.9, 88.9);
        test_exp2::<f32>(-150.0, 128.0);
        test_exp10::<f32>(-45.9, 38.9);
    }

    #[test]
    fn test_f64() {
        test_exp::<f64>(-745.9, 709.9);
        test_exp_m1::<f64>(-708.9, 709.9);
        test_exp2::<f64>(-1075.0, 1024.0);
        test_exp10::<f64>(-323.9, 308.9);
    }
}
