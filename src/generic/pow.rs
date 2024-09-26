use super::exp::{exp_inner_common, exp_split};
use super::log::{hi_lo_log_special_poly, log_split};
use super::{int_is_odd, is_int, is_odd_int, Exp, Log};
use crate::double::{DenormDouble, SemiDouble};
use crate::traits::{CastInto as _, Int as _};

pub(crate) fn pow<F: Log + Exp>(x: F, y: F) -> F {
    let (nx, xedelta) = x.normalize_arg();
    let (ny, _) = y.normalize_arg();
    let xexp = nx.raw_exp();
    let yexp = ny.raw_exp();

    if yexp == F::RawExp::ZERO || nx == F::one() {
        // pow(x, 0) = 1
        // pow(1, y) = 1
        F::one()
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
            if ny.sign() {
                // pow(±0, y) = ±inf when y < 0
                F::INFINITY.copysign(nx)
            } else {
                // pow(±0, y) = ±0 when y > 0
                nx
            }
        } else {
            // y is not an odd integer
            if ny.sign() {
                // pow(±0, y) = +inf when y < 0
                F::INFINITY
            } else {
                // pow(±0, y) = 0 when y > 0
                F::ZERO
            }
        }
    } else if xexp == F::MAX_RAW_EXP {
        // x = ±inf
        if nx.sign() {
            // x = -inf
            if is_odd_int(ny) {
                // y is an odd integer
                if ny.sign() {
                    // pow(-inf, y) = -0 when y < 0
                    -F::ZERO
                } else {
                    // pow(-inf, y) = -inf when y > 0
                    F::neg_infinity()
                }
            } else {
                // y is not an odd integer
                if ny.sign() {
                    // pow(-inf, y) = -0 when y < 0
                    F::ZERO
                } else {
                    // pow(-inf, y) = inf when y > 0
                    F::INFINITY
                }
            }
        } else {
            // x = +inf
            if ny.sign() {
                // pow(+inf, y) = 0 when y < 0
                F::ZERO
            } else {
                // pow(+inf, y) = +inf when y > 0
                F::INFINITY
            }
        }
    } else if yexp == F::MAX_RAW_EXP {
        // y = ±inf
        if nx == -F::one() {
            // pow(-1, ±inf) = 1
            F::one()
        } else if ny.sign() {
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
    } else if nx.sign() && !is_int(ny) {
        // pow(x, y) = NaN when x < 0 and y is finite and not integer
        F::NAN
    } else {
        // logx = log(|x|)
        let logx = hi_lo_log_inner(nx.abs(), xedelta).to_semi();

        // ylx = y * log(|x|)
        let ylx = (logx * y).to_norm();

        // |z| = |x|^y = exp(y * log(|x|))
        let absz = if ylx.hi() >= F::exp_hi_th() {
            F::INFINITY
        } else if ylx.hi() <= F::exp_lo_th() {
            F::ZERO
        } else {
            let (k, r_hi, r_lo) = exp_split(ylx.hi());
            let r_lo = r_lo + ylx.lo();

            exp_inner_common(k, r_hi, r_lo)
        };

        if nx.sign() && int_is_odd(ny) {
            -absz
        } else {
            absz
        }
    }
}

pub(super) fn hi_lo_log_inner<F: Log>(x: F, edelta: F::Exp) -> DenormDouble<F> {
    // Algorithm based on one used by the msun math library:
    //  * log(1 + r) = p * s + 2 * s
    //  * s = r / (2 + r)
    //  * p = (log(1 + s) - log(1 - s) - 2 * s) / s

    // Split x * 2^edelta = 2^k * (1 + r)
    //  - k is an integer
    //  - sqrt(2) / 2 <= 1 + r < sqrt(2)
    let (k, r) = log_split(x, edelta);

    // rp2 = 2 + r
    let rp2 = SemiDouble::new_qadd11(F::two(), r);

    // s = r / (2 + r)
    let s = (SemiDouble::new(r) / rp2).to_semi();
    let s2 = s.square().to_semi();

    // p = (log(1 + s) - log(1 - s) - 2 * s) / s
    let p = hi_lo_log_special_poly(s2).to_semi();

    // t1 = k * log(2)
    let kf: F = k.cast_into();
    let t1 = DenormDouble::new(F::ln_2_hi(), F::ln_2_lo()).pmul1(kf);

    // t2 = log(1 + r) = p * s + 2 * s
    let ps = p * s;
    let twos = s.pmul1(F::two());
    let t2 = twos.to_denorm().qadd2(ps);

    // log(2^k * (1 + r)) = t1 + t2
    t1.qadd2(t2)
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test<F: Float + FloatMath>() {
        use crate::pow;

        let f = F::parse;

        assert_is_nan!(pow(F::NAN, F::NAN));
        assert_is_nan!(pow(F::ZERO, F::NAN));
        assert_is_nan!(pow(-F::ZERO, F::NAN));
        assert_is_nan!(pow(F::two(), F::NAN));
        assert_is_nan!(pow(F::INFINITY, F::NAN));
        assert_is_nan!(pow(F::neg_infinity(), F::NAN));
        assert_is_nan!(pow(F::NAN, F::one()));
        assert_is_nan!(pow(F::NAN, F::INFINITY));
        assert_is_nan!(pow(F::NAN, F::neg_infinity()));
        assert_is_nan!(pow(f("-3"), f("0.5")));
        assert_total_eq!(pow(F::ZERO, f("-33")), F::INFINITY);
        assert_total_eq!(pow(-F::ZERO, f("-33")), F::neg_infinity());
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
        assert_total_eq!(pow(F::ZERO, F::neg_infinity()), F::INFINITY);
        assert_total_eq!(pow(-F::ZERO, F::neg_infinity()), F::INFINITY);
        assert_total_eq!(pow(F::one(), F::ZERO), F::one());
        assert_total_eq!(pow(F::one(), -F::ZERO), F::one());
        assert_total_eq!(pow(F::one(), f("33")), F::one());
        assert_total_eq!(pow(F::one(), f("-33")), F::one());
        assert_total_eq!(pow(F::one(), f("33.5")), F::one());
        assert_total_eq!(pow(F::one(), f("-33.5")), F::one());
        assert_total_eq!(pow(F::one(), f("34")), F::one());
        assert_total_eq!(pow(F::one(), f("-34.0")), F::one());
        assert_total_eq!(pow(F::one(), F::INFINITY), F::one());
        assert_total_eq!(pow(F::one(), F::neg_infinity()), F::one());
        assert_total_eq!(pow(F::one(), F::NAN), F::one());
        assert_total_eq!(pow(-F::one(), F::INFINITY), F::one());
        assert_total_eq!(pow(-F::one(), F::neg_infinity()), F::one());
        assert_total_eq!(pow(f("0.5"), F::INFINITY), F::ZERO);
        assert_total_eq!(pow(f("0.5"), F::neg_infinity()), F::INFINITY);
        assert_total_eq!(pow(f("-0.5"), F::INFINITY), F::ZERO);
        assert_total_eq!(pow(f("-0.5"), F::neg_infinity()), F::INFINITY);
        assert_total_eq!(pow(f("1.5"), F::INFINITY), F::INFINITY);
        assert_total_eq!(pow(f("1.5"), F::neg_infinity()), F::ZERO);
        assert_total_eq!(pow(f("-1.5"), F::INFINITY), F::INFINITY);
        assert_total_eq!(pow(f("-1.5"), F::neg_infinity()), F::ZERO);
        assert_total_eq!(pow(F::INFINITY, F::ZERO), F::one());
        assert_total_eq!(pow(F::INFINITY, -F::ZERO), F::one());
        assert_total_eq!(pow(F::INFINITY, f("33")), F::INFINITY);
        assert_total_eq!(pow(F::INFINITY, f("-33")), F::ZERO);
        assert_total_eq!(pow(F::INFINITY, f("33.5")), F::INFINITY);
        assert_total_eq!(pow(F::INFINITY, f("-33.5")), F::ZERO);
        assert_total_eq!(pow(F::INFINITY, f("34.0")), F::INFINITY);
        assert_total_eq!(pow(F::INFINITY, f("-34.0")), F::ZERO);
        assert_total_eq!(pow(F::neg_infinity(), F::ZERO), F::one());
        assert_total_eq!(pow(F::neg_infinity(), -F::ZERO), F::one());
        assert_total_eq!(pow(F::neg_infinity(), f("33")), F::neg_infinity());
        assert_total_eq!(pow(F::neg_infinity(), f("-33")), -F::ZERO);
        assert_total_eq!(pow(F::neg_infinity(), f("33.5")), F::INFINITY);
        assert_total_eq!(pow(F::neg_infinity(), f("-33.5")), F::ZERO);
        assert_total_eq!(pow(F::neg_infinity(), f("34.0")), F::INFINITY);
        assert_total_eq!(pow(F::neg_infinity(), f("-34.0")), F::ZERO);
        assert_total_eq!(pow(F::two(), F::two()), f("4"));
        assert_total_eq!(pow(F::two(), -F::two()), f("0.25"));
        assert_total_eq!(pow(-F::two(), f("3")), f("-8"));
        assert_total_eq!(pow(-F::two(), f("-3")), f("-0.125"));
        assert_total_eq!(pow(f("3.5"), f("3")), f("42.875"));
        assert_total_eq!(pow(f("10"), f("4")), f("10000"));
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
