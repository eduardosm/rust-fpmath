use super::exp::exp_inner_common;
use super::{Exp, round_as_i_f};
use crate::traits::Int as _;

pub(crate) trait Exp2: Exp {
    const LN_2: Self;
    const EXP2_LO_TH: Self;
    const EXP2_HI_TH: Self;
}

/// Returns 2 raised to `x`.
pub(crate) fn exp2<F: Exp2>(x: F) -> F {
    if x >= F::EXP2_HI_TH {
        // also handles x = inf
        F::INFINITY
    } else if x <= F::EXP2_LO_TH {
        // also handles x = -inf
        F::ZERO
    } else {
        let e = x.raw_exp();
        if e == F::RawExp::ZERO {
            // x is zero or subnormal
            // 2^x ~= 1
            F::ONE
        } else if e == F::MAX_RAW_EXP {
            // x is NaN, propagate
            x
        } else {
            exp2_inner(x)
        }
    }
}

/// Calculates `2^x` where:
///
///  * `x` is not zero, subnormal, nan nor infinity
///  * `x` is less than `EXP2_HI_TH`
fn exp2_inner<F: Exp2>(x: F) -> F {
    // Split x into k, r_hi, r_lo such as:
    //  - x = k + (r_hi + r_lo)*log2(e)
    //  - k is an integer
    //  - |r_hi| <= 0.5*ln(2)
    let (k, r_hi, r_lo) = exp2_split(x);

    // Calculate 2^x = exp(k*ln(2) + r_hi + r_lo)
    exp_inner_common(k, r_hi, r_lo)
}

/// Splits `x` into `(k, r_hi, r_lo)`
///
/// Such as:
/// * `x = k + (r_hi + r_lo)*log2(e)`
/// * `k` is an integer
/// * `|r_hi| <= 0.5*ln(2)`
#[inline]
fn exp2_split<F: Exp2>(x: F) -> (i32, F, F) {
    let (k, kf) = round_as_i_f(x);
    let t = x - kf;

    let (t_hi, t_lo) = t.split_hi_lo();
    let r_hi = t_hi * F::LN_2_HI;
    let r_lo = t_hi * F::LN_2_LO + t_lo * F::LN_2;

    (k, r_hi, r_lo)
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::Float;

    fn test<F: Float + FloatMath>(lo_th: &str, hi_th: &str) {
        use crate::exp2;

        let f = F::parse;

        let lo_th = f(lo_th);
        let hi_th = f(hi_th);

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

    #[test]
    fn test_f32() {
        test::<f32>("-150", "128");
    }

    #[test]
    fn test_f64() {
        test::<f64>("-1075", "1024");
    }
}
