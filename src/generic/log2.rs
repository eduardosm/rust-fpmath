use super::{log::log_split, Log};
use crate::traits::{CastInto as _, Int as _, Like};

pub(crate) trait Log2<L = Like<Self>>: Log {
    fn log2_e_hi() -> Self;
    fn log2_e_lo() -> Self;
}

/// Calculates `log2(x)`
///
/// `x` must be normal and positive.
fn log2_inner<F: Log2>(x: F, edelta: F::Exp) -> F {
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
    let hr2 = F::half() * r * r;
    let t1_hi = (r - hr2).split_hi();
    let t1_lo = ((r - t1_hi) - hr2) + s * (hr2 + p);

    // t2 = log2(1 + r) = log(1 + r) * log2(e) = t1 * log2(e)
    let t2_hi = t1_hi * F::log2_e_hi();
    let t2_lo = (t1_hi + t1_lo) * F::log2_e_lo() + t1_lo * F::log2_e_hi();

    // log2(x) = k + log2(1 + r) = k + t2
    let kf: F = k.cast_into();
    let t4_hi = (kf + t2_hi).purify();
    let t4_lo = t2_lo + ((kf - t4_hi) + t2_hi);
    t4_hi + t4_lo
}

pub(crate) fn log2<F: Log2>(x: F) -> F {
    let (y, edelta) = x.normalize_arg();
    let yexp = y.raw_exp();
    if yexp == F::RawExp::ZERO {
        // log2(±0) = -inf
        F::neg_infinity()
    } else if y.sign() {
        // x < 0, log2(x) = NaN
        F::NAN
    } else if yexp == F::MAX_RAW_EXP {
        if y.raw_mant() == F::Raw::ZERO {
            // log2(inf) = inf
            F::INFINITY
        } else {
            // NaN, propagate
            y
        }
    } else {
        log2_inner(y, edelta)
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test<F: Float + FloatMath>() {
        use crate::log2;

        assert_is_nan!(log2(F::NAN));
        assert_is_nan!(log2(-F::one()));
        assert_is_nan!(log2(F::neg_infinity()));
        assert_total_eq!(log2(F::ZERO), F::neg_infinity());
        assert_total_eq!(log2(-F::ZERO), F::neg_infinity());
        assert_total_eq!(log2(F::INFINITY), F::INFINITY);
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
