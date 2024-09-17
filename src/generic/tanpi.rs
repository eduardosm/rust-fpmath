use super::{reduce_half_mul_pi, tan::tan_inner, ReduceHalfMulPi, Tan};
use crate::traits::{CastFrom as _, Int as _};

pub(crate) fn tanpi<F: ReduceHalfMulPi + Tan>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        // tanpi(inf or NaN) = NaN
        F::NAN
    } else if e == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // tanpi(±0) = ±0
        x
    } else if e <= F::RawExp::from(F::MANT_BITS) {
        // very small: tanpi(x) ~= x * π

        // scale temporarily to avoid temporary subnormal numbers
        let logscale = F::Exp::TWO * F::Exp::cast_from(F::MANT_BITS);
        let scale = F::exp2i_fast(logscale);
        let descale = F::exp2i_fast(-logscale);

        let (x_hi, x_lo) = (x * scale).split_hi_lo();
        let y_hi = x_hi * F::pi_hi();
        let y_lo = x_hi * F::pi_lo() + x_lo * F::pi();
        (y_hi + y_lo) * descale
    } else {
        let (n, y_hi, y_lo) = reduce_half_mul_pi(x);
        let inv = (n & 1) != 0;
        if inv && y_hi == F::ZERO {
            F::INFINITY.copysign(x)
        } else {
            tan_inner(y_hi, y_lo, inv)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test<F: Float + FloatMath>() {
        use crate::tanpi;

        assert_is_nan!(tanpi(F::NAN));
        assert_is_nan!(tanpi(F::INFINITY));
        assert_is_nan!(tanpi(F::neg_infinity()));
        assert_total_eq!(tanpi(F::ZERO), F::ZERO);
        assert_total_eq!(tanpi(-F::ZERO), -F::ZERO);
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
