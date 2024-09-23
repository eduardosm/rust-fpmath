use super::exp::exp_inner_common;
use super::{round_as_i_f, Exp};
use crate::traits::{Int as _, Like};

pub(crate) trait Exp2<L = Like<Self>>: Exp {
    fn ln_2() -> Self;
    fn exp2_lo_th() -> Self;
    fn exp2_hi_th() -> Self;
}

/// Returns 2 raised to `x`.
pub(crate) fn exp2<F: Exp2>(x: F) -> F {
    if x >= F::exp2_hi_th() {
        // also handles x = inf
        F::INFINITY
    } else if x <= F::exp2_lo_th() {
        // also handles x = -inf
        F::ZERO
    } else {
        let e = x.raw_exp();
        if e == F::RawExp::ZERO {
            // x is zero or subnormal
            // 2^x ~= 1
            F::one()
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
    let r_hi = t_hi * F::ln_2_hi();
    let r_lo = t_hi * F::ln_2_lo() + t_lo * F::ln_2();

    (k, r_hi, r_lo)
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test<F: Float + FloatMath>(lo_th: &str, hi_th: &str) {
        use crate::exp2;

        let f = F::parse;

        let lo_th = f(lo_th);
        let hi_th = f(hi_th);

        assert_is_nan!(exp2(F::NAN));
        assert_total_eq!(exp2(F::INFINITY), F::INFINITY);
        assert_total_eq!(exp2(F::neg_infinity()), F::ZERO);
        assert_total_eq!(exp2(F::ZERO), F::one());
        assert_total_eq!(exp2(-F::ZERO), F::one());
        assert_total_eq!(exp2(F::one()), F::two());
        assert_total_eq!(exp2(F::two()), f("4"));
        assert_total_eq!(exp2(f("32")), f("4294967296"));
        assert_total_eq!(exp2(-F::one()), F::half());
        assert_total_eq!(exp2(-F::two()), f("0.25"));
        assert_total_eq!(exp2(f("-3")), f("0.125"));
        assert_total_eq!(exp2(f("-4")), f("0.0625"));
        assert_total_eq!(exp2(lo_th), F::ZERO);
        assert_total_eq!(exp2(lo_th - F::one()), F::ZERO);
        assert_total_eq!(exp2(lo_th - F::two()), F::ZERO);
        assert_total_eq!(exp2(hi_th), F::INFINITY);
        assert_total_eq!(exp2(hi_th + F::one()), F::INFINITY);
        assert_total_eq!(exp2(hi_th + F::two()), F::INFINITY);
    }

    #[test]
    fn test_f32() {
        test::<f32>("-150", "128");
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f32() {
        test::<crate::SoftF32>("-150", "128");
    }

    #[test]
    fn test_f64() {
        test::<f64>("-1075", "1024");
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f64() {
        test::<crate::SoftF64>("-1075", "1024");
    }
}
