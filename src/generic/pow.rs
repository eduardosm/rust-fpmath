use super::{int_is_odd, is_int, is_odd_int};
use crate::traits::{Float, Int as _};

pub(crate) trait Pow: Float {
    fn pow_finite(x: Self, xedelta: Self::Exp, y: Self, sign: bool) -> Self;

    fn powi_finite(x: Self, xedelta: Self::Exp, y: i32) -> Self;
}

pub(crate) fn pow<F: Pow>(x: F, y: F) -> F {
    let (nx, xedelta) = x.normalize_arg();
    let (ny, _) = y.normalize_arg();
    let xexp = nx.raw_exp();
    let yexp = ny.raw_exp();

    if yexp == F::RawExp::ZERO || nx == F::ONE {
        // pow(x, 0) = 1
        // pow(1, y) = 1
        F::ONE
    } else if (yexp == F::MAX_RAW_EXP && ny.raw_mant() != F::Raw::ZERO)
        || (xexp == F::MAX_RAW_EXP && nx.raw_mant() != F::Raw::ZERO)
    {
        // pow(x, NaN) = NaN when x != 1
        // pow(NaN, y) = NaN when y != 0
        F::NAN
    } else if xexp == F::RawExp::ZERO {
        // x = ±0
        if is_odd_int(ny) {
            // y is an odd integer
            if ny.is_sign_negative() {
                // pow(±0, y) = ±inf when y < 0
                F::INFINITY.copysign(nx)
            } else {
                // pow(±0, y) = ±0 when y > 0
                nx
            }
        } else {
            // y is not an odd integer
            if ny.is_sign_negative() {
                // pow(±0, y) = +inf when y < 0
                F::INFINITY
            } else {
                // pow(±0, y) = 0 when y > 0
                F::ZERO
            }
        }
    } else if xexp == F::MAX_RAW_EXP {
        // x = ±inf
        if nx.is_sign_negative() {
            // x = -inf
            if is_odd_int(ny) {
                // y is an odd integer
                if ny.is_sign_negative() {
                    // pow(-inf, y) = -0 when y < 0
                    -F::ZERO
                } else {
                    // pow(-inf, y) = -inf when y > 0
                    F::NEG_INFINITY
                }
            } else {
                // y is not an odd integer
                if ny.is_sign_negative() {
                    // pow(-inf, y) = -0 when y < 0
                    F::ZERO
                } else {
                    // pow(-inf, y) = inf when y > 0
                    F::INFINITY
                }
            }
        } else {
            // x = +inf
            if ny.is_sign_negative() {
                // pow(+inf, y) = 0 when y < 0
                F::ZERO
            } else {
                // pow(+inf, y) = +inf when y > 0
                F::INFINITY
            }
        }
    } else if yexp == F::MAX_RAW_EXP {
        // y = ±inf
        if nx == -F::ONE {
            // pow(-1, ±inf) = 1
            F::ONE
        } else if ny.is_sign_negative() {
            // y = -inf
            if xexp < F::EXP_OFFSET {
                // pow(x, -inf) = inf when |x| < 1
                F::INFINITY
            } else {
                // pow(x, -inf) = 0 when |x| > 1
                F::ZERO
            }
        } else {
            // y = +inf
            if xexp < F::EXP_OFFSET {
                // pow(x, +inf) = 0 when |x| < 1
                F::ZERO
            } else {
                // pow(x, +inf) = inf when |x| > 1
                F::INFINITY
            }
        }
    } else if nx.is_sign_negative() && !is_int(ny) {
        // pow(x, y) = NaN when x < 0 and y is finite and not integer
        F::NAN
    } else {
        F::pow_finite(nx, xedelta, y, x.is_sign_negative() && int_is_odd(y))
    }
}

pub(crate) fn powi<F: Pow>(x: F, y: i32) -> F {
    let (nx, xedelta) = x.normalize_arg();
    let xexp = nx.raw_exp();

    if y == 0 || nx == F::ONE {
        // pow(x, 0) = 1
        // pow(1, y) = 1
        F::ONE
    } else if xexp == F::MAX_RAW_EXP && nx.raw_mant() != F::Raw::ZERO {
        // pow(NaN, y) = NaN when y != 0
        F::NAN
    } else if xexp == F::RawExp::ZERO {
        // x = ±0
        if (y & 1) != 0 {
            // y is odd
            if y < 0 {
                // pow(±0, y) = ±inf when y < 0
                F::INFINITY.copysign(nx)
            } else {
                // pow(±0, y) = ±0 when y > 0
                nx
            }
        } else {
            // y is even
            if y < 0 {
                // pow(±0, y) = +inf when y < 0
                F::INFINITY
            } else {
                // pow(±0, y) = 0 when y > 0
                F::ZERO
            }
        }
    } else if xexp == F::MAX_RAW_EXP {
        // x = ±inf
        if nx.is_sign_negative() {
            // x = -inf
            if (y & 1) != 0 {
                // y is odd
                if y < 0 {
                    // pow(-inf, y) = -0 when y < 0
                    -F::ZERO
                } else {
                    // pow(-inf, y) = -inf when y > 0
                    F::NEG_INFINITY
                }
            } else {
                // y is even
                if y < 0 {
                    // pow(-inf, y) = -0 when y < 0
                    F::ZERO
                } else {
                    // pow(-inf, y) = inf when y > 0
                    F::INFINITY
                }
            }
        } else {
            // x = +inf
            if y < 0 {
                // pow(+inf, y) = 0 when y < 0
                F::ZERO
            } else {
                // pow(+inf, y) = +inf when y > 0
                F::INFINITY
            }
        }
    } else {
        F::powi_finite(nx, xedelta, y)
    }
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::Float;

    fn test_pow<F: Float + FloatMath>() {
        use crate::pow;

        let f = F::parse;

        assert_is_nan!(pow(F::NAN, F::NAN));
        assert_is_nan!(pow(F::ZERO, F::NAN));
        assert_is_nan!(pow(-F::ZERO, F::NAN));
        assert_is_nan!(pow(F::TWO, F::NAN));
        assert_is_nan!(pow(F::INFINITY, F::NAN));
        assert_is_nan!(pow(F::NEG_INFINITY, F::NAN));
        assert_is_nan!(pow(F::NAN, F::ONE));
        assert_is_nan!(pow(F::NAN, F::INFINITY));
        assert_is_nan!(pow(F::NAN, F::NEG_INFINITY));
        assert_is_nan!(pow(f("-3"), f("0.5")));
        assert_total_eq!(pow(F::ZERO, f("-33")), F::INFINITY);
        assert_total_eq!(pow(-F::ZERO, f("-33")), F::NEG_INFINITY);
        assert_total_eq!(pow(F::ZERO, f("-33.5")), F::INFINITY);
        assert_total_eq!(pow(-F::ZERO, f("-33.5")), F::INFINITY);
        assert_total_eq!(pow(F::ZERO, f("-34")), F::INFINITY);
        assert_total_eq!(pow(-F::ZERO, f("-34")), F::INFINITY);
        assert_total_eq!(pow(F::ZERO, f("33")), F::ZERO);
        assert_total_eq!(pow(-F::ZERO, f("33")), -F::ZERO);
        assert_total_eq!(pow(F::ZERO, f("33.5")), F::ZERO);
        assert_total_eq!(pow(-F::ZERO, f("33.5")), F::ZERO);
        assert_total_eq!(pow(F::ZERO, f("34.0")), F::ZERO);
        assert_total_eq!(pow(-F::ZERO, f("34.0")), F::ZERO);
        assert_total_eq!(pow(F::ZERO, F::INFINITY), F::ZERO);
        assert_total_eq!(pow(-F::ZERO, F::INFINITY), F::ZERO);
        assert_total_eq!(pow(F::ZERO, F::NEG_INFINITY), F::INFINITY);
        assert_total_eq!(pow(-F::ZERO, F::NEG_INFINITY), F::INFINITY);
        assert_total_eq!(pow(F::ONE, F::ZERO), F::ONE);
        assert_total_eq!(pow(F::ONE, -F::ZERO), F::ONE);
        assert_total_eq!(pow(F::ONE, f("33")), F::ONE);
        assert_total_eq!(pow(F::ONE, f("-33")), F::ONE);
        assert_total_eq!(pow(F::ONE, f("33.5")), F::ONE);
        assert_total_eq!(pow(F::ONE, f("-33.5")), F::ONE);
        assert_total_eq!(pow(F::ONE, f("34")), F::ONE);
        assert_total_eq!(pow(F::ONE, f("-34.0")), F::ONE);
        assert_total_eq!(pow(F::ONE, F::INFINITY), F::ONE);
        assert_total_eq!(pow(F::ONE, F::NEG_INFINITY), F::ONE);
        assert_total_eq!(pow(F::ONE, F::NAN), F::ONE);
        assert_total_eq!(pow(-F::ONE, F::INFINITY), F::ONE);
        assert_total_eq!(pow(-F::ONE, F::NEG_INFINITY), F::ONE);
        assert_total_eq!(pow(f("0.5"), F::INFINITY), F::ZERO);
        assert_total_eq!(pow(f("0.5"), F::NEG_INFINITY), F::INFINITY);
        assert_total_eq!(pow(f("-0.5"), F::INFINITY), F::ZERO);
        assert_total_eq!(pow(f("-0.5"), F::NEG_INFINITY), F::INFINITY);
        assert_total_eq!(pow(f("1.5"), F::INFINITY), F::INFINITY);
        assert_total_eq!(pow(f("1.5"), F::NEG_INFINITY), F::ZERO);
        assert_total_eq!(pow(f("-1.5"), F::INFINITY), F::INFINITY);
        assert_total_eq!(pow(f("-1.5"), F::NEG_INFINITY), F::ZERO);
        assert_total_eq!(pow(F::INFINITY, F::ZERO), F::ONE);
        assert_total_eq!(pow(F::INFINITY, -F::ZERO), F::ONE);
        assert_total_eq!(pow(F::INFINITY, f("33")), F::INFINITY);
        assert_total_eq!(pow(F::INFINITY, f("-33")), F::ZERO);
        assert_total_eq!(pow(F::INFINITY, f("33.5")), F::INFINITY);
        assert_total_eq!(pow(F::INFINITY, f("-33.5")), F::ZERO);
        assert_total_eq!(pow(F::INFINITY, f("34.0")), F::INFINITY);
        assert_total_eq!(pow(F::INFINITY, f("-34.0")), F::ZERO);
        assert_total_eq!(pow(F::NEG_INFINITY, F::ZERO), F::ONE);
        assert_total_eq!(pow(F::NEG_INFINITY, -F::ZERO), F::ONE);
        assert_total_eq!(pow(F::NEG_INFINITY, f("33")), F::NEG_INFINITY);
        assert_total_eq!(pow(F::NEG_INFINITY, f("-33")), -F::ZERO);
        assert_total_eq!(pow(F::NEG_INFINITY, f("33.5")), F::INFINITY);
        assert_total_eq!(pow(F::NEG_INFINITY, f("-33.5")), F::ZERO);
        assert_total_eq!(pow(F::NEG_INFINITY, f("34.0")), F::INFINITY);
        assert_total_eq!(pow(F::NEG_INFINITY, f("-34.0")), F::ZERO);
        assert_total_eq!(pow(F::TWO, F::TWO), f("4"));
        assert_total_eq!(pow(F::TWO, -F::TWO), f("0.25"));
        assert_total_eq!(pow(-F::TWO, f("3")), f("-8"));
        assert_total_eq!(pow(-F::TWO, f("-3")), f("-0.125"));
        assert_total_eq!(pow(f("3.5"), f("3")), f("42.875"));
        assert_total_eq!(pow(f("10"), f("4")), f("10000"));
    }

    fn test_powi<F: Float + FloatMath>() {
        use crate::powi;

        let f = F::parse;

        assert_is_nan!(powi(F::NAN, 1));
        assert_total_eq!(powi(F::ZERO, -33), F::INFINITY);
        assert_total_eq!(powi(-F::ZERO, -33), F::NEG_INFINITY);
        assert_total_eq!(powi(F::ZERO, -34), F::INFINITY);
        assert_total_eq!(powi(-F::ZERO, -34), F::INFINITY);
        assert_total_eq!(powi(F::ZERO, 33), F::ZERO);
        assert_total_eq!(powi(-F::ZERO, 33), -F::ZERO);
        assert_total_eq!(powi(F::ZERO, 34), F::ZERO);
        assert_total_eq!(powi(-F::ZERO, 34), F::ZERO);
        assert_total_eq!(powi(F::ONE, 0), F::ONE);
        assert_total_eq!(powi(F::ONE, 33), F::ONE);
        assert_total_eq!(powi(F::ONE, -33), F::ONE);
        assert_total_eq!(powi(F::ONE, 34), F::ONE);
        assert_total_eq!(powi(F::ONE, -34), F::ONE);
        assert_total_eq!(powi(F::INFINITY, 0), F::ONE);
        assert_total_eq!(powi(F::INFINITY, 33), F::INFINITY);
        assert_total_eq!(powi(F::INFINITY, -33), F::ZERO);
        assert_total_eq!(powi(F::INFINITY, 34), F::INFINITY);
        assert_total_eq!(powi(F::INFINITY, -34), F::ZERO);
        assert_total_eq!(powi(F::NEG_INFINITY, 0), F::ONE);
        assert_total_eq!(powi(F::NEG_INFINITY, -0), F::ONE);
        assert_total_eq!(powi(F::NEG_INFINITY, 33), F::NEG_INFINITY);
        assert_total_eq!(powi(F::NEG_INFINITY, -33), -F::ZERO);
        assert_total_eq!(powi(F::NEG_INFINITY, 34), F::INFINITY);
        assert_total_eq!(powi(F::NEG_INFINITY, -34), F::ZERO);
        assert_total_eq!(powi(F::TWO, 2), f("4"));
        assert_total_eq!(powi(F::TWO, -2), f("0.25"));
        assert_total_eq!(powi(-F::TWO, 3), f("-8"));
        assert_total_eq!(powi(-F::TWO, -3), f("-0.125"));
        assert_total_eq!(powi(f("3.5"), 3), f("42.875"));
        assert_total_eq!(powi(f("10"), 4), f("10000"));
    }

    #[test]
    fn test_f32() {
        test_pow::<f32>();
        test_powi::<f32>();
    }

    #[test]
    fn test_f64() {
        test_pow::<f64>();
        test_powi::<f64>();
    }
}
