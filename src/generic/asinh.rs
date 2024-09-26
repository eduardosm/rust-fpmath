use super::log::{log_hi_lo_inner, log_inner};
use super::sqrt::hi_lo_sqrt_hi_lo_inner;
use super::Log;
use crate::double::DenormDouble;
use crate::traits::Int as _;

pub(crate) fn asinh<F: Log>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP || e <= (F::EXP_OFFSET - F::RawExp::from(F::MANT_BITS)) {
        // asinh(±inf) = ±inf
        // or
        // propagate NaN
        // or
        // very small, includes subnormal and zero
        // asinh(x) ~= x
        // also handles asinh(-0) = -0
        x
    } else if e > (F::RawExp::from(F::MANT_BITS) + F::EXP_OFFSET) {
        let y = log_inner(x.abs(), F::Exp::ONE);
        y.copysign(x)
    } else {
        asinh_inner(x)
    }
}

fn asinh_inner<F: Log>(x: F) -> F {
    let absx = x.abs();
    let x2 = x * x;

    // t1 = x^2 + 1
    let t1 = DenormDouble::new_add11(x2, F::one());

    // t2 = sqrt(x^2 + 1)
    let t2 = hi_lo_sqrt_hi_lo_inner(t1);

    // t3 = |x| + sqrt(x^2 + 1)
    let t3 = t2.qadd1(absx);

    // t4 = |asinh(x)| = log(|x| + sqrt(x^2 + 1))
    let t4 = log_hi_lo_inner(t3.hi(), t3.lo());

    // asinh(x) = sgn(x) * |asinh(x)|
    t4.copysign(x)
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test<F: Float + FloatMath>() {
        use crate::asinh;

        assert_is_nan!(asinh(F::NAN));
        assert_total_eq!(asinh(F::INFINITY), F::INFINITY);
        assert_total_eq!(asinh(F::neg_infinity()), F::neg_infinity());
        assert_total_eq!(asinh(F::ZERO), F::ZERO);
        assert_total_eq!(asinh(-F::ZERO), -F::ZERO);
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
