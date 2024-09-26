use super::log::{log_hi_lo_inner, log_inner};
use super::sqrt::hi_lo_sqrt_hi_lo_inner;
use super::Log;
use crate::double::DenormDouble;
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
    let t1 = if x < F::two() {
        // y = x - 1
        let y = x - F::one();
        let y2 = y * y;
        let twoy = F::two() * y;

        // t1 = x^2 - 1 = y^2 + 2 * y
        DenormDouble::new_qadd11(twoy, y2)
    } else {
        let x2 = x * x;

        // t1 = x^2 - 1
        DenormDouble::new_qsub11(x2, F::one())
    };

    // t2 = sqrt(x^2 - 1)
    let t2 = hi_lo_sqrt_hi_lo_inner(t1);

    // t3 = x + sqrt(x^2 - 1)
    let t3 = t2.qradd1(x);

    // acosh(x) = log(x + sqrt(x^2 - 1))
    log_hi_lo_inner(t3.hi(), t3.lo())
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
