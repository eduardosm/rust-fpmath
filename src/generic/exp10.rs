use super::exp::exp_inner_common;
use super::{Exp, round_as_i_f};
use crate::traits::Int as _;

pub(crate) trait Exp10: Exp {
    const LOG2_10: Self;
    const LOG10_2_HI: Self;
    const LOG10_2_LO: Self;
    const LN_10: Self;
    const LN_10_HI: Self;
    const LN_10_LO: Self;
    const EXP10_LO_TH: Self;
    const EXP10_HI_TH: Self;
}

/// Returns 10 raised to `x`.
pub(crate) fn exp10<F: Exp10>(x: F) -> F {
    if x >= F::EXP10_HI_TH {
        // also handles x = inf
        F::INFINITY
    } else if x <= F::EXP10_LO_TH {
        // also handles x = -inf
        F::ZERO
    } else {
        let e = x.raw_exp();
        if e == F::RawExp::ZERO {
            // x is zero or subnormal
            // 10^x ~= 1
            F::ONE
        } else if e == F::MAX_RAW_EXP {
            // x is NaN, propagate
            x
        } else {
            exp10_inner(x)
        }
    }
}

fn exp10_inner<F: Exp10>(x: F) -> F {
    // Split x into k, r_hi, r_lo such as:
    //  - x = k*log10(2) + (r_hi + r_lo)*log10(e)
    //  - k is an integer
    //  - |r| <= 0.5*ln(2)
    let (k, r_hi, r_lo) = exp10_split(x);

    // Calculate 10^x = exp(k*ln(2) + r_hi + r_lo)
    exp_inner_common(k, r_hi, r_lo)
}

/// Splits `x` into `(k, r_hi, r_lo)`
///
/// Such as:
/// * `x = k*log10(2) + (r_hi + r_lo)*log10(e)`
/// * `k` is an integer
/// * `|r| <= 0.5*ln(2)`
#[inline]
fn exp10_split<F: Exp10>(x: F) -> (i32, F, F) {
    let y = x * F::LOG2_10;
    let (k, kf) = round_as_i_f(y);
    // `kf * LOG10_2_HI` is exact because the lower bits
    // of `LOG10_2_HI` are zero.
    let t_hi = x - kf * F::LOG10_2_HI;
    let t_lo = -kf * F::LOG10_2_LO;
    let (t_hi, t_lo) = F::norm_hi_lo_splitted(t_hi, t_lo);

    let r_hi = t_hi * F::LN_10_HI;
    let r_lo = t_hi * F::LN_10_LO + t_lo * F::LN_10;

    (k, r_hi, r_lo)
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::Float;

    fn test<F: Float + FloatMath>(lo_th: &str, hi_th: &str) {
        use crate::exp10;

        let f = F::parse;

        let lo_th = f(lo_th);
        let hi_th = f(hi_th);

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
        test::<f32>("-45.9", "38.9");
    }

    #[test]
    fn test_f64() {
        test::<f64>("-323.9", "308.9");
    }
}
