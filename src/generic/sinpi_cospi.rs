use super::sin_cos::{cos_inner, sin_inner};
use super::{reduce_half_mul_pi, ReduceHalfMulPi, SinCos};
use crate::traits::{CastFrom as _, Int as _};

pub(crate) fn sinpi<F: SinCos + ReduceHalfMulPi>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        // sinpi(inf or nan) = nan
        F::NAN
    } else if e == F::RawExp::ZERO {
        if x.raw_mant() == F::Raw::ZERO {
            // sinpi(±0) = ±0
            x
        } else {
            // subnormal: sinpi(x) ~= x * π

            // scale temporarily to avoid temporary subnormal numbers
            let logscale = F::Exp::TWO * F::Exp::cast_from(F::MANT_BITS);
            let scale = F::exp2i_fast(logscale);
            let descale = F::exp2i_fast(-logscale);

            let (x_hi, x_lo) = (x * scale).split_hi_lo();
            let y_hi = x_hi * F::pi_hi();
            let y_lo = x_hi * F::pi_lo() + x_lo * F::pi();
            (y_hi + y_lo) * descale
        }
    } else {
        let (n, y_hi, y_lo) = reduce_half_mul_pi(x);

        match n {
            0 => sin_inner(y_hi, y_lo),
            1 => cos_inner(y_hi, y_lo),
            2 => -sin_inner(y_hi, y_lo),
            3 => -cos_inner(y_hi, y_lo),
            _ => unreachable!(),
        }
    }
}

pub(crate) fn cospi<F: SinCos + ReduceHalfMulPi>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        // cospi(inf or nan) = nan
        F::NAN
    } else if e <= F::RawExp::from(F::MANT_BITS) {
        // subnormal or zero, cospi(x) ~= 1
        F::one()
    } else {
        let (n, y_hi, y_lo) = reduce_half_mul_pi(x);

        match n {
            0 => cos_inner(y_hi, y_lo),
            1 => -sin_inner(y_hi, y_lo),
            2 => -cos_inner(y_hi, y_lo),
            3 => sin_inner(y_hi, y_lo),
            _ => unreachable!(),
        }
    }
}

pub(crate) fn sinpi_cospi<F: SinCos + ReduceHalfMulPi>(x: F) -> (F, F) {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        // sinpi(inf or nan) = nan
        // cospi(inf or nan) = nan
        (F::NAN, F::NAN)
    } else if e == F::RawExp::ZERO {
        if x.raw_mant() == F::Raw::ZERO {
            // sinpi(±0) = ±0
            // cospi(±0) = 1
            (x, F::one())
        } else {
            // subnormal:
            // sinpi(x) ~= x * π
            // cospi(x) ~= 1

            // scale temporarily to avoid temporary subnormal numbers
            let logscale = F::Exp::TWO * F::Exp::cast_from(F::MANT_BITS);
            let scale = F::exp2i_fast(logscale);
            let descale = F::exp2i_fast(-logscale);

            let (x_hi, x_lo) = (x * scale).split_hi_lo();
            let y_hi = x_hi * F::pi_hi();
            let y_lo = x_hi * F::pi_lo() + x_lo * F::pi();
            ((y_hi + y_lo) * descale, F::one())
        }
    } else {
        let (n, y_hi, y_lo) = reduce_half_mul_pi(x);

        let sin = sin_inner(y_hi, y_lo);
        let cos = cos_inner(y_hi, y_lo);
        match n {
            0 => (sin, cos),
            1 => (cos, -sin),
            2 => (-sin, -cos),
            3 => (-cos, sin),
            _ => unreachable!(),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test<F: Float + FloatMath>() {
        use crate::{cospi, sinpi, sinpi_cospi};

        let test_nan = |arg: F| {
            let sin1 = sinpi(arg);
            let cos1 = cospi(arg);
            let (sin2, cos2) = sinpi_cospi(arg);
            assert_is_nan!(sin1);
            assert_is_nan!(cos1);
            assert_is_nan!(sin2);
            assert_is_nan!(cos2);
        };

        let test_value = |arg: F, expected_sin: F, expected_cos: F| {
            let sin1 = sinpi(arg);
            let cos1 = cospi(arg);
            let (sin2, cos2) = sinpi_cospi(arg);
            assert_total_eq!(sin1, expected_sin);
            assert_total_eq!(cos1, expected_cos);
            assert_total_eq!(sin2, expected_sin);
            assert_total_eq!(cos2, expected_cos);
        };

        test_nan(F::NAN);
        test_nan(F::INFINITY);
        test_nan(F::neg_infinity());
        test_value(F::ZERO, F::ZERO, F::one());
        test_value(-F::ZERO, -F::ZERO, F::one());
    }

    #[test]
    fn test_f32() {
        test::<f32>();
    }

    #[test]
    fn test_soft_f32() {
        test::<crate::SoftF32>();
    }

    #[test]
    fn test_f64() {
        test::<f64>();
    }

    #[test]
    fn test_soft_f64() {
        test::<crate::SoftF64>();
    }
}
