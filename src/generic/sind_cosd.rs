use super::sin_cos::{cos_inner, sin_inner};
use super::{reduce_90_deg, Reduce90Deg, SinCos};

pub(crate) fn sind<F: SinCos + Reduce90Deg>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        // sind(inf or nan) = nan
        F::NAN
    } else if e <= F::RawExp::from(F::MANT_BITS) {
        // subnormal or zero, sind(x) ~= x * (π/180)
        // also handles sind(-0) = -0
        x * F::deg_to_rad()
    } else {
        let (n, y_hi, y_lo) = reduce_90_deg(x);

        match n {
            0 => sin_inner(y_hi, y_lo),
            1 => cos_inner(y_hi, y_lo),
            2 => -sin_inner(y_hi, y_lo),
            3 => -cos_inner(y_hi, y_lo),
            _ => unreachable!(),
        }
    }
}

pub(crate) fn cosd<F: SinCos + Reduce90Deg>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        // cosd(inf or nan) = nan
        F::NAN
    } else if e <= F::RawExp::from(F::MANT_BITS) {
        // subnormal or zero, cosd(x) ~= 1
        F::one()
    } else {
        let (n, y_hi, y_lo) = reduce_90_deg(x);

        match n {
            0 => cos_inner(y_hi, y_lo),
            1 => -sin_inner(y_hi, y_lo),
            2 => -cos_inner(y_hi, y_lo),
            3 => sin_inner(y_hi, y_lo),
            _ => unreachable!(),
        }
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
        (x * F::deg_to_rad(), F::one())
    } else {
        let (n, y_hi, y_lo) = reduce_90_deg(x);

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
        use crate::{cosd, sind, sind_cosd};

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
