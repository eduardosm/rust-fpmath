use super::sin_cos::{cos_quadrant, sin_cos_quadrant, sin_quadrant};
use super::{ReduceHalfMulPi, SinCos, reduce_half_mul_pi};
use crate::double::SemiDouble;
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

            let sx = SemiDouble::new(x * scale);
            let y = sx * F::PI_EX;
            y.to_single() * descale
        }
    } else {
        let (n, y) = reduce_half_mul_pi(x);
        sin_quadrant(n, x.sign(), y.hi(), y.lo())
    }
}

pub(crate) fn cospi<F: SinCos + ReduceHalfMulPi>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        // cospi(inf or nan) = nan
        F::NAN
    } else if e <= F::RawExp::from(F::MANT_BITS) {
        // subnormal or zero, cospi(x) ~= 1
        F::ONE
    } else {
        let (n, y) = reduce_half_mul_pi(x);
        cos_quadrant(n, y.hi(), y.lo())
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
            (x, F::ONE)
        } else {
            // subnormal:
            // sinpi(x) ~= x * π
            // cospi(x) ~= 1

            // scale temporarily to avoid temporary subnormal numbers
            let logscale = F::Exp::TWO * F::Exp::cast_from(F::MANT_BITS);
            let scale = F::exp2i_fast(logscale);
            let descale = F::exp2i_fast(-logscale);

            let sx = SemiDouble::new(x * scale);
            let y = sx * F::PI_EX;
            (y.to_single() * descale, F::ONE)
        }
    } else {
        let (n, y) = reduce_half_mul_pi(x);
        sin_cos_quadrant(n, x.sign(), y.hi(), y.lo())
    }
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::Float;

    fn test<F: Float + FloatMath>() {
        use crate::{cospi, sinpi, sinpi_cospi};

        let f = F::parse;

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
        test_nan(F::NEG_INFINITY);
        test_value(F::ZERO, F::ZERO, F::ONE);
        test_value(-F::ZERO, -F::ZERO, F::ONE);
        test_value(f("0.5"), F::ONE, F::ZERO);
        test_value(f("-0.5"), -F::ONE, F::ZERO);
        test_value(f("1.0"), F::ZERO, -F::ONE);
        test_value(f("-1.0"), -F::ZERO, -F::ONE);
        test_value(f("1.5"), -F::ONE, F::ZERO);
        test_value(f("-1.5"), F::ONE, F::ZERO);
        test_value(f("2.0"), F::ZERO, F::ONE);
        test_value(f("-2.0"), -F::ZERO, F::ONE);
        test_value(f("2.5"), F::ONE, F::ZERO);
        test_value(f("-2.5"), -F::ONE, F::ZERO);
    }

    #[test]
    fn test_f32() {
        test::<f32>();
    }

    #[test]
    fn test_f64() {
        test::<f64>();
    }
}
