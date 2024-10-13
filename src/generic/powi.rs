use super::exp::exp_split;
use super::log::hi_lo_log_inner;
use super::{scalbn, Exp, Log};
use crate::double::{DenormDouble, SemiDouble};
use crate::traits::{CastInto as _, Int as _};

pub(crate) fn powi<F: Log + Exp>(x: F, y: i32) -> F {
    let (nx, xedelta) = x.normalize_arg();
    let xexp = nx.raw_exp();

    if y == 0 || nx == F::one() {
        // pow(x, 0) = 1
        // pow(1, y) = 1
        F::one()
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
        if nx.sign() {
            // x = -inf
            if (y & 1) != 0 {
                // y is odd
                if y < 0 {
                    // pow(-inf, y) = -0 when y < 0
                    -F::ZERO
                } else {
                    // pow(-inf, y) = -inf when y > 0
                    F::neg_infinity()
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
        // logx = log(|x|)
        let logx = hi_lo_log_inner(nx.abs(), xedelta).to_semi();

        // Split y = sum(y_i)
        // |x|^y = prod(|x|^y_i)
        //       = prod(e^(log(|x|)*y_i))
        //       = prod(e^(r_i)) * 2^(sum(k_i))
        // log(|x|) * y_i = r_i + k_i*ln(2)

        // k_total = sum(k_i)
        // z = prod(e^(r_i))
        let mut k_total = 0;
        let mut z = SemiDouble::one();

        let absy = y.unsigned_abs();
        let mut yshift = 0;
        while yshift < 32 {
            // yf is a chunk of y (one of the y_i above)
            let mask: u32 = ((F::MANT_MASK << 1) | F::Raw::ONE).cast_into();
            let yf: F = ((absy >> yshift) & mask).cast_into();
            let yf = (yf * F::exp2i_fast(yshift.cast_into())).set_sign(y < 0);
            yshift += F::MANT_BITS + 1;

            // ylx = yf * log(|x|)
            let ylx = (logx * yf).to_norm();

            if ylx.hi() >= F::exp_hi_th() {
                if (absy & 1) == 0 || !nx.sign() {
                    return F::INFINITY;
                } else {
                    return F::neg_infinity();
                }
            } else if ylx.hi() <= F::exp_lo_th() {
                if (absy & 1) == 0 || !nx.sign() {
                    return F::ZERO;
                } else {
                    return -F::ZERO;
                }
            } else {
                let (k, r_hi, r_lo) = exp_split(ylx.hi());

                // t = |x|^yf / 2^k = exp(yf * log(|x|)) / 2^k
                let r = DenormDouble::new(r_hi, r_lo + ylx.lo());
                let t = hi_lo_exp_inner(r).to_semi();

                k_total += k;
                z = (z * t).to_semi();
            }
        }

        let absz = scalbn(z.to_single(), k_total);
        if nx.sign() && (y & 1) != 0 {
            -absz
        } else {
            absz
        }
    }
}

fn hi_lo_exp_inner<F: Exp>(r: DenormDouble<F>) -> DenormDouble<F> {
    // Calculates exp(r_hi + r_lo)
    // Similar to `exp_inner_common` in exp.rs, but returns hi/lo parts
    // and assumes k=0

    let r_single = r.to_single();
    let r2 = r_single * r_single;

    // t1 = 2 - 2 * r / (exp(r) - 1)
    let t1 = r_single + F::exp_special_poly(r2);

    // t2 = r * t1 / (2 - t1)
    let t2 = r_single * t1 / (F::two() - t1);

    // t3 = r + t2
    let t3 = r.qadd1(t2);

    // 1 + t3
    t3.qradd1(F::one())
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test<F: Float + FloatMath>() {
        use crate::powi;

        let f = F::parse;

        assert_is_nan!(powi(F::NAN, 1));
        assert_total_eq!(powi(F::ZERO, -33), F::INFINITY);
        assert_total_eq!(powi(-F::ZERO, -33), F::neg_infinity());
        assert_total_eq!(powi(F::ZERO, -34), F::INFINITY);
        assert_total_eq!(powi(-F::ZERO, -34), F::INFINITY);
        assert_total_eq!(powi(F::ZERO, 33), F::ZERO);
        assert_total_eq!(powi(-F::ZERO, 33), -F::ZERO);
        assert_total_eq!(powi(F::ZERO, 34), F::ZERO);
        assert_total_eq!(powi(-F::ZERO, 34), F::ZERO);
        assert_total_eq!(powi(F::one(), 0), F::one());
        assert_total_eq!(powi(F::one(), 33), F::one());
        assert_total_eq!(powi(F::one(), -33), F::one());
        assert_total_eq!(powi(F::one(), 34), F::one());
        assert_total_eq!(powi(F::one(), -34), F::one());
        assert_total_eq!(powi(F::INFINITY, 0), F::one());
        assert_total_eq!(powi(F::INFINITY, 33), F::INFINITY);
        assert_total_eq!(powi(F::INFINITY, -33), F::ZERO);
        assert_total_eq!(powi(F::INFINITY, 34), F::INFINITY);
        assert_total_eq!(powi(F::INFINITY, -34), F::ZERO);
        assert_total_eq!(powi(F::neg_infinity(), 0), F::one());
        assert_total_eq!(powi(F::neg_infinity(), -0), F::one());
        assert_total_eq!(powi(F::neg_infinity(), 33), F::neg_infinity());
        assert_total_eq!(powi(F::neg_infinity(), -33), -F::ZERO);
        assert_total_eq!(powi(F::neg_infinity(), 34), F::INFINITY);
        assert_total_eq!(powi(F::neg_infinity(), -34), F::ZERO);
        assert_total_eq!(powi(F::two(), 2), f("4"));
        assert_total_eq!(powi(F::two(), -2), f("0.25"));
        assert_total_eq!(powi(-F::two(), 3), f("-8"));
        assert_total_eq!(powi(-F::two(), -3), f("-0.125"));
        assert_total_eq!(powi(f("3.5"), 3), f("42.875"));
        assert_total_eq!(powi(f("10"), 4), f("10000"));
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
