use super::log::{log_hi_lo_inner, log_inner};
use super::sqrt::hi_lo_sqrt_hi_lo_inner;
use super::Log;
use crate::traits::Int as _;

pub(crate) fn acosh<F: Log>(x: F) -> F {
    let e = x.raw_exp();
    if x < F::one() {
        // x < 1, acosh(x) is NaN
        F::NAN
    } else if x == F::one() {
        // asinh(1) = 0
        F::ZERO
    } else if e == F::MAX_RAW_EXP {
        // x is infinity or NaN
        // acosh(x) = x
        x
    } else if e > (F::RawExp::from(F::MANT_BITS) + F::EXP_OFFSET) {
        log_inner(x, F::Exp::ONE)
    } else {
        acosh_inner(x)
    }
}

fn acosh_inner<F: Log>(x: F) -> F {
    // t1 = x^2 - 1
    let (t1_hi, t1_lo) = if x < F::two() {
        // y = x - 1
        let y = x - F::one();
        let y2 = y * y;
        let twoy = F::two() * y;

        // t1 = x^2 - 1 = y^2 + 2 * y
        let t1_hi = (y2 + twoy).purify();
        let t1_lo = (twoy - t1_hi) + y2;

        (t1_hi, t1_lo)
    } else {
        let x2 = x * x;

        // t1 = x^2 - 1
        let t1_hi = (x2 - F::one()).purify();
        let t1_lo = (x2 - t1_hi) - F::one();

        (t1_hi, t1_lo)
    };

    // t2 = sqrt(x^2 - 1)
    let (t2_hi, t2_lo) = hi_lo_sqrt_hi_lo_inner(t1_hi, t1_lo);

    // t3 = x + sqrt(x^2 - 1)
    let t3_hi = (x + t2_hi).purify();
    let t3_lo = ((x - t3_hi) + t2_hi) + t2_lo;

    // acosh(x) = log(x + sqrt(x^2 - 1))
    log_hi_lo_inner(t3_hi, t3_lo)
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test<F: Float + FloatMath>() {
        use crate::acosh;

        assert_is_nan!(acosh(F::NAN));
        assert_is_nan!(acosh(F::neg_infinity()));
        assert_is_nan!(acosh(-F::one()));
        assert_is_nan!(acosh(F::ZERO));
        assert_is_nan!(acosh(F::half()));
        assert_total_eq!(acosh(F::INFINITY), F::INFINITY);
        assert_total_eq!(acosh(F::one()), F::ZERO);
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
