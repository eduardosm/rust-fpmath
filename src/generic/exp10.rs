use super::exp::exp_inner_common;
use super::{round_as_i_f, Exp};
use crate::traits::{Int as _, Like};

pub(crate) trait Exp10<L = Like<Self>>: Exp {
    fn log2_10() -> Self;
    fn log10_2_hi() -> Self;
    fn log10_2_lo() -> Self;
    fn ln_10() -> Self;
    fn ln_10_hi() -> Self;
    fn ln_10_lo() -> Self;
    fn exp10_lo_th() -> Self;
    fn exp10_hi_th() -> Self;
}

/// Splits `x` into `(k, r_hi, r_lo)`
///
/// Such as:
/// * `x = k*log10(2) + (r_hi + r_lo)*log10(e)`
/// * `k` is an integer
/// * `|r| <= 0.5*ln(2)`
#[inline]
fn exp10_split<F: Exp10>(x: F) -> (i32, F, F) {
    let y = x * F::log2_10();
    let (k, kf) = round_as_i_f(y);
    // `kf * LOG10_2_HI` is exact because the lower bits
    // of `LOG10_2_HI` are zero.
    let t_hi = x - kf * F::log10_2_hi();
    let t_lo = -kf * F::log10_2_lo();
    let (t_hi, t_lo) = F::norm_hi_lo_splitted(t_hi, t_lo);

    let r_hi = t_hi * F::ln_10_hi();
    let r_lo = t_hi * F::ln_10_lo() + t_lo * F::ln_10();

    (k, r_hi, r_lo)
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

/// Returns 10 raised to `x`.
pub(crate) fn exp10<F: Exp10>(x: F) -> F {
    if x >= F::exp10_hi_th() {
        // also handles x = inf
        F::INFINITY
    } else if x <= F::exp10_lo_th() {
        // also handles x = -inf
        F::ZERO
    } else {
        let e = x.raw_exp();
        if e == F::RawExp::ZERO {
            // x is zero or subnormal
            // 10^x ~= 1
            F::one()
        } else if e == F::MAX_RAW_EXP {
            // x is NaN, propagate
            x
        } else {
            exp10_inner(x)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test<F: Float + FloatMath>(lo_th: &str, hi_th: &str) {
        use crate::exp10;

        let f = F::parse;

        let lo_th = f(lo_th);
        let hi_th = f(hi_th);

        assert_is_nan!(exp10(F::NAN));
        assert_total_eq!(exp10(F::INFINITY), F::INFINITY);
        assert_total_eq!(exp10(F::neg_infinity()), F::ZERO);
        assert_total_eq!(exp10(F::ZERO), F::one());
        assert_total_eq!(exp10(-F::ZERO), F::one());
        assert_total_eq!(exp10(F::one()), f("10"));
        assert_total_eq!(exp10(F::two()), f("100"));
        assert_total_eq!(exp10(lo_th), F::ZERO);
        assert_total_eq!(exp10(lo_th - F::one()), F::ZERO);
        assert_total_eq!(exp10(lo_th - F::two()), F::ZERO);
        assert_total_eq!(exp10(hi_th), F::INFINITY);
        assert_total_eq!(exp10(hi_th + F::one()), F::INFINITY);
        assert_total_eq!(exp10(hi_th + F::two()), F::INFINITY);
    }

    #[test]
    fn test_f32() {
        test::<f32>("-45.9", "38.9");
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f32() {
        test::<crate::SoftF32>("-45.9", "38.9");
    }

    #[test]
    fn test_f64() {
        test::<f64>("-323.9", "308.9");
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f64() {
        test::<crate::SoftF64>("-323.9", "308.9");
    }
}
