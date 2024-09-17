use super::exp::{exp_inner_common, exp_split};
use super::log::{hi_lo_log_special_poly, log_split};
use super::{Exp, Log};
use crate::traits::{CastInto as _, Float, Int as _};

pub(super) fn hi_lo_log_inner<F: Log>(x: F, edelta: F::Exp) -> (F, F) {
    // Algorithm based on one used by the msun math library:
    //  * log(1 + r) = p * s + 2 * s
    //  * s = r / (2 + r)
    //  * p = (log(1 + s) - log(1 - s) - 2 * s) / s

    // Split x * 2^edelta = 2^k * (1 + r)
    //  - k is an integer
    //  - sqrt(2) / 2 <= 1 + r < sqrt(2)
    let (k, r) = log_split(x, edelta);
    let (r_hi, r_lo) = r.split_hi_lo();

    // rp2 = 2 + r
    let rp2_hi = (F::two() + r).split_hi();
    let rp2_lo = (F::two() - rp2_hi) + r;

    // s = r / (2 + r)
    let (s_hi, s_lo) = F::div_hi_lo(r_hi, r_lo, rp2_hi, rp2_lo);
    let (s_hi, s_lo) = F::norm_hi_lo_splitted(s_hi, s_lo);

    let s2_hi = s_hi * s_hi;
    let s2_lo = F::two() * s_hi * s_lo + s_lo * s_lo;
    let (s2_hi, s2_lo) = F::norm_hi_lo_splitted(s2_hi, s2_lo);

    // p = (log(1 + s) - log(1 - s) - 2 * s) / s
    let (p_hi, p_lo) = hi_lo_log_special_poly(s2_hi, s2_lo);
    let (p_hi, p_lo) = F::norm_hi_lo_splitted(p_hi, p_lo);

    // t1 = k * log(2)
    let kf: F = k.cast_into();
    let t1_hi = kf * F::ln_2_hi();
    let t1_lo = kf * F::ln_2_lo();

    // t2 = log(1 + r) = p * s + 2 * s
    let ps_hi = p_hi * s_hi;
    let ps_lo = p_hi * s_lo + p_lo * s_hi + p_lo * s_lo;

    let twos_hi = F::two() * s_hi;
    let twos_lo = F::two() * s_lo;

    let t2_hi = (ps_hi + twos_hi).purify();
    let t2_lo = (((twos_hi - t2_hi) + ps_hi) + twos_lo) + ps_lo;

    // t3 = log(2^k * (1 + r)) = t1 + t2
    let t3_hi = (t1_hi + t2_hi).purify();
    let t3_lo = (((t1_hi - t3_hi) + t2_hi) + t1_lo) + t2_lo;

    (t3_hi, t3_lo)
}

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
        let (logx_hi, logx_lo) = hi_lo_log_inner(nx.abs(), xedelta);
        let (logx_hi, logx_lo) = F::norm_hi_lo_splitted(logx_hi, logx_lo);

        let (y_hi, y_lo) = y.split_hi_lo();

        // ylx = y * log(|x|)
        let ylx_hi = logx_hi * y_hi;
        let ylx_lo = logx_hi * y_lo + logx_lo * y;
        let (ylx_hi, ylx_lo) = F::norm_hi_lo_full(ylx_hi, ylx_lo);

        // |z| = |x|^y = exp(y * log(|x|))
        let absz = if ylx_hi >= F::exp_hi_th() {
            F::INFINITY
        } else if ylx_hi <= F::exp_lo_th() {
            F::ZERO
        } else {
            let (k, r_hi, r_lo) = exp_split(ylx_hi);
            let r_lo = r_lo + ylx_lo;

            exp_inner_common(k, r_hi, r_lo)
        };

        if nx.sign() && is_odd(ny) {
            -absz
        } else {
            absz
        }
    }
}

fn is_int<F: Float>(x: F) -> bool {
    let e = x.raw_exp();
    if e > F::EXP_OFFSET + F::RawExp::from(F::MANT_BITS) {
        true
    } else if e < F::EXP_OFFSET {
        false
    } else {
        let frac_shift = (F::EXP_OFFSET + F::RawExp::from(F::MANT_BITS)) - e;
        (x.to_raw() & !(F::Raw::MAX << frac_shift)) == F::Raw::ZERO
    }
}

fn is_odd_int<F: Float>(x: F) -> bool {
    let e = x.raw_exp();
    if e > F::EXP_OFFSET + F::RawExp::from(F::MANT_BITS) || e < F::EXP_OFFSET {
        // infinity, an even integer or only fractional part (less than 1)
        false
    } else {
        let frac_shift = (F::EXP_OFFSET + F::RawExp::from(F::MANT_BITS)) - e;
        if (x.to_raw() & !(F::Raw::MAX << frac_shift)) != F::Raw::ZERO {
            // not an integer
            false
        } else {
            ((x.mant() >> frac_shift) & F::Raw::ONE) == F::Raw::ONE
        }
    }
}

// like `is_odd_int`, but assumes that `x` is an integer
fn is_odd<F: Float>(x: F) -> bool {
    let e = x.raw_exp();
    if e > F::EXP_OFFSET + F::RawExp::from(F::MANT_BITS) {
        false
    } else {
        let frac_shift = (F::EXP_OFFSET + F::RawExp::from(F::MANT_BITS)) - e;
        ((x.mant() >> frac_shift) & F::Raw::ONE) == F::Raw::ONE
    }
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

    #[test]
    fn test_soft_f32() {
        test::<crate::SoftF32>();
    }

    #[test]
    fn test_f64() {
        test::<f64>();
    }

    #[test]
    fn test_soft_f64() {
        test::<crate::SoftF64>();
    }
}
