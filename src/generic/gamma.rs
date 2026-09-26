use super::is_int;
use crate::traits::{Float, Int as _};

pub(crate) trait Gamma: Float {
    fn gamma_finite(x: Self) -> Self;

    fn ln_gamma_finite(x: Self) -> (Self, i8);
}

pub(crate) fn gamma<F: Gamma>(x: F) -> F {
    let e = x.raw_exp();
    let sign = x.is_sign_negative();
    if e == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // gamma(±0) = ±inf
        F::INFINITY.copysign(x)
    } else if e == F::MAX_RAW_EXP {
        if !sign && x.raw_mant() == F::Raw::ZERO {
            // gamma(inf) = inf
            F::INFINITY
        } else {
            // gamma(NaN or -inf) = NaN
            F::NAN
        }
    } else if sign && is_int(x) {
        // gamma(neg integer) = NaN
        F::NAN
    } else {
        F::gamma_finite(x)
    }
}

pub(crate) fn ln_gamma<F: Gamma>(x: F) -> (F, i8) {
    let e = x.raw_exp();
    let sign = x.is_sign_negative();
    if e == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // ln_gamma(0) = inf
        (F::INFINITY, if sign { -1 } else { 1 })
    } else if e == F::MAX_RAW_EXP {
        if !sign && x.raw_mant() == F::Raw::ZERO {
            // ln_gamma(inf) = inf
            (F::INFINITY, 1)
        } else {
            // ln_gamma(NaN or -inf) = NaN
            (F::NAN, 0)
        }
    } else if sign && is_int(x) {
        // ln_gamma(neg integer) = inf
        (F::INFINITY, 0)
    } else if x == F::ONE || x == F::TWO {
        // ln_gamma(1 or 2) = 0
        // ensure positive zero
        (F::ZERO, 1)
    } else {
        F::ln_gamma_finite(x)
    }
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::Float;

    fn test_gamma<F: Float + FloatMath>(underflow_firsts: &[F]) {
        use crate::gamma;

        assert_is_nan!(gamma(F::NAN));
        assert_is_nan!(gamma(F::NEG_INFINITY));
        assert_is_nan!(gamma(-F::ONE));
        assert_is_nan!(gamma(-F::TWO));
        assert_is_nan!(gamma(-F::LARGEST));
        assert_total_eq!(gamma(F::INFINITY), F::INFINITY);
        assert_total_eq!(gamma(F::ZERO), F::INFINITY);
        assert_total_eq!(gamma(-F::ZERO), F::NEG_INFINITY);
        assert_total_eq!(gamma(F::ONE), F::ONE);
        assert_total_eq!(gamma(F::TWO), F::ONE);

        // Once the result underflows to zero it must still carry the sign of
        // the true (tiny, non-zero) value, which alternates between
        // consecutive negative integers.
        for &first in underflow_firsts {
            let mut x = first;
            let mut negative = true;
            for _ in 0..8 {
                assert_total_eq!(gamma(x), if negative { -F::ZERO } else { F::ZERO });
                x = x - F::ONE;
                negative = !negative;
            }
        }
    }

    fn test_ln_gamma<F: Float + FloatMath>() {
        use crate::ln_gamma;

        let test_nan = |x: F| {
            let (r, sign) = ln_gamma(x);
            assert_is_nan!(r);
            assert_eq!(sign, 0);
        };
        let test_value = |x: F, r: F, sign: i8| {
            let (res, res_sign) = ln_gamma(x);
            assert_total_eq!(res, r);
            assert_eq!(res_sign, sign);
        };

        test_nan(F::NAN);
        test_nan(F::NEG_INFINITY);
        test_value(F::INFINITY, F::INFINITY, 1);
        test_value(F::ZERO, F::INFINITY, 1);
        test_value(-F::ZERO, F::INFINITY, -1);
        test_value(-F::ONE, F::INFINITY, 0);
        test_value(-F::TWO, F::INFINITY, 0);
        test_value(-F::LARGEST, F::INFINITY, 0);
        test_value(F::ONE, F::ZERO, 1);
        test_value(F::TWO, F::ZERO, 1);
    }

    #[test]
    fn test_f32() {
        test_gamma::<f32>(&[-42.5, -42.25, -42.75]);
        test_ln_gamma::<f32>();
    }

    #[test]
    fn test_f64() {
        test_gamma::<f64>(&[-180.5, -180.25, -180.75]);
        test_ln_gamma::<f64>();
    }
}
