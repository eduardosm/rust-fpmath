use super::sin_cos::{cos_quadrant, sin_cos_quadrant, sin_quadrant};
use super::{Reduce90Deg, SinCos, reduce_90_deg};

pub(crate) fn sind<F: SinCos + Reduce90Deg>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        // sind(inf or nan) = nan
        F::NAN
    } else if e <= F::RawExp::from(F::MANT_BITS) {
        // subnormal or zero, sind(x) ~= x * (π/180)
        // also handles sind(-0) = -0
        x * F::DEG_TO_RAD
    } else {
        let (n, y) = reduce_90_deg(x);
        sin_quadrant(n, x.sign(), y.hi(), y.lo())
    }
}

pub(crate) fn cosd<F: SinCos + Reduce90Deg>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        // cosd(inf or nan) = nan
        F::NAN
    } else if e <= F::RawExp::from(F::MANT_BITS) {
        // subnormal or zero, cosd(x) ~= 1
        F::ONE
    } else {
        let (n, y) = reduce_90_deg(x);
        cos_quadrant(n, y.hi(), y.lo())
    }
}

pub(crate) fn sind_cosd<F: SinCos + Reduce90Deg>(x: F) -> (F, F) {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        // sind(inf or nan) = nan
        // cosd(inf or nan) = nan
        (F::NAN, F::NAN)
    } else if e <= F::RawExp::from(F::MANT_BITS) {
        // subnormal or zero
        // sind(x) ~= x * (π/180)
        // cosd(x) ~= 1
        // also handles sind(-0) = -0
        (x * F::DEG_TO_RAD, F::ONE)
    } else {
        let (n, y) = reduce_90_deg(x);
        sin_cos_quadrant(n, x.sign(), y.hi(), y.lo())
    }
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::Float;

    fn test<F: Float + FloatMath>() {
        use crate::{cosd, sind, sind_cosd};

        let f = F::parse;

        let test_nan = |arg: F| {
            let sin1 = sind(arg);
            let cos1 = cosd(arg);
            let (sin2, cos2) = sind_cosd(arg);
            assert_is_nan!(sin1);
            assert_is_nan!(cos1);
            assert_is_nan!(sin2);
            assert_is_nan!(cos2);
        };

        let test_value = |arg: F, expected_sin: F, expected_cos: F| {
            let sin1 = sind(arg);
            let cos1 = cosd(arg);
            let (sin2, cos2) = sind_cosd(arg);
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
        test_value(f("90"), F::ONE, F::ZERO);
        test_value(f("-90"), -F::ONE, F::ZERO);
        test_value(f("180"), F::ZERO, -F::ONE);
        test_value(f("-180"), -F::ZERO, -F::ONE);
        test_value(f("270"), -F::ONE, F::ZERO);
        test_value(f("-270"), F::ONE, F::ZERO);
        test_value(f("360"), F::ZERO, F::ONE);
        test_value(f("-360"), -F::ZERO, F::ONE);
        test_value(f("450"), F::ONE, F::ZERO);
        test_value(f("-450"), -F::ONE, F::ZERO);
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
