use super::{log::log_split, Log};
use crate::double::{DenormDouble, SemiDouble};
use crate::traits::{CastInto as _, Int as _, Like};

pub(crate) trait Log10<L = Like<Self>>: Log {
    fn log10_e_ex() -> SemiDouble<Self>;
    fn log10_2_hi() -> Self;
    fn log10_2_lo() -> Self;
}

pub(crate) fn log10<F: Log10>(x: F) -> F {
    let (y, edelta) = x.normalize_arg();
    let yexp = y.raw_exp();
    if yexp == F::RawExp::ZERO {
        // log10(±0) = -inf
        F::neg_infinity()
    } else if y.sign() {
        // x < 0, log10(x) = NaN
        F::NAN
    } else if yexp == F::MAX_RAW_EXP {
        if y.raw_mant() == F::Raw::ZERO {
            // log10(inf) = inf
            F::INFINITY
        } else {
            // NaN, propagate
            y
        }
    } else {
        log10_inner(y, edelta)
    }
}

/// Calculates `log10(x)`
///
/// `x` must be finite normal and positive.
fn log10_inner<F: Log10>(x: F, edelta: F::Exp) -> F {
    // Algorithm based on one used by the msun math library:
    //  * log(1 + r) = p * s + 2 * s
    //  * s = r / (2 + r)
    //  * p = (log(1 + s) - log(1 - s) - 2 * s) / s

    // Split x * 2^edelta = 2^k * (1 + r)
    //  - k is an integer
    //  - sqrt(2) / 2 <= 1 + r < sqrt(2)
    let (k, r) = log_split(x, edelta);

    // s = r / (2 + r)
    // So, log(1 + r) = log(1 + s) - log(1 - s)
    let s = r / (F::two() + r);

    // p = (log(1 + s) - log(1 - s) - 2 * s) / s
    let p = F::log_special_poly(s);

    // t1 = log(1 + r) = p * s + 2 * s
    //    = r - s * (r - p)
    //    = r - (0.5 * r^2 - s * (0.5 * r^2 + p))
    // Split t1 into t1_hi + t1_lo for better accuracy
    let hr2 = (F::half() * r * r).purify();
    let t1 = DenormDouble::new_qsub11(r, hr2)
        .qadd1(s * (hr2 + p))
        .to_semi();

    // t2 = log10(1 + r) = log(1 + r) * log10(e) = t1 * log10(e)
    let t2 = t1 * F::log10_e_ex();

    // t3 = k * log10(2)
    let kf: F = k.cast_into();
    let t3 = DenormDouble::new(F::log10_2_hi(), F::log10_2_lo()).pmul1(kf);

    // log10(x) = k * log10(2) + log10(1 + r) = t2 + t3
    let t4 = t3.qadd2(t2);
    t4.to_single()
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test<F: Float + FloatMath>() {
        use crate::log10;

        assert_is_nan!(log10(F::NAN));
        assert_is_nan!(log10(-F::one()));
        assert_is_nan!(log10(F::neg_infinity()));
        assert_total_eq!(log10(F::ZERO), F::neg_infinity());
        assert_total_eq!(log10(-F::ZERO), F::neg_infinity());
        assert_total_eq!(log10(F::INFINITY), F::INFINITY);
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
