use super::asin_acos::{acos_inner, asin_inner};
use super::{AsinAcos, DivPi};
use crate::double::SemiDouble;
use crate::traits::{CastFrom as _, Int as _};

pub(crate) fn asinpi<F: AsinAcos + DivPi>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::EXP_OFFSET && x.raw_mant() == F::Raw::ZERO {
        // asinpi(±1) = ±0.5
        F::HALF.copysign(x)
    } else if e >= F::EXP_OFFSET {
        // NaN or |x| > 1 (including infinity)
        F::NAN
    } else if e == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // asinpi(±0) = ±0
        x
    } else if e <= F::RawExp::from(F::MANT_BITS) {
        // very small, asinpi(x) ~= x / π

        // scale temporarily to avoid temporary subnormal numbers
        let logscale = F::Exp::TWO * F::Exp::cast_from(F::MANT_BITS);
        let scale = F::exp2i_fast(logscale);
        let descale = F::exp2i_fast(-logscale);

        let nx = SemiDouble::new(x * scale);

        (nx * F::FRAC_1_PI_EX).to_single() * descale
    } else {
        let y = asin_inner(x).to_semi();

        (y * F::FRAC_1_PI_EX).to_single()
    }
}

pub(crate) fn acospi<F: AsinAcos + DivPi>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::EXP_OFFSET && x.raw_mant() == F::Raw::ZERO {
        if x.sign() {
            // acospi(-1) = 1
            F::ONE
        } else {
            // acospi(1) = 0
            F::ZERO
        }
    } else if e >= F::EXP_OFFSET {
        // NaN or |x| > 1 (including infinity)
        F::NAN
    } else if e == F::RawExp::ZERO {
        // subnormal or zero
        // acospi(x) ~= 0.5
        F::HALF
    } else {
        let y = acos_inner(x).to_semi();

        (y * F::FRAC_1_PI_EX).to_single()
    }
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::Float;

    fn test_asinpi<F: Float + FloatMath>() {
        use crate::asinpi;

        let f = F::parse;

        assert_is_nan!(asinpi(F::NAN));
        assert_is_nan!(asinpi(f("1.5")));
        assert_is_nan!(asinpi(f("-1.5")));
        assert_is_nan!(asinpi(F::INFINITY));
        assert_is_nan!(asinpi(F::NEG_INFINITY));
        assert_total_eq!(asinpi(F::ZERO), F::ZERO);
        assert_total_eq!(asinpi(-F::ZERO), -F::ZERO);
        assert_total_eq!(asinpi(F::ONE), F::HALF);
        assert_total_eq!(asinpi(-F::ONE), -F::HALF);
    }

    fn test_acospi<F: Float + FloatMath>() {
        use crate::acospi;

        let f = F::parse;

        assert_is_nan!(acospi(F::NAN));
        assert_is_nan!(acospi(f("1.5")));
        assert_is_nan!(acospi(f("-1.5")));
        assert_is_nan!(acospi(F::INFINITY));
        assert_is_nan!(acospi(F::NEG_INFINITY));
        assert_total_eq!(acospi(F::ZERO), F::HALF);
        assert_total_eq!(acospi(-F::ZERO), F::HALF);
        assert_total_eq!(acospi(F::ONE), F::ZERO);
        assert_total_eq!(acospi(-F::ONE), F::ONE);
    }

    #[test]
    fn test_f32() {
        test_asinpi::<f32>();
        test_acospi::<f32>();
    }

    #[test]
    fn test_f64() {
        test_asinpi::<f64>();
        test_acospi::<f64>();
    }
}
