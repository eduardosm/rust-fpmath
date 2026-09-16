use super::Ln;
use super::ln::ln_hi_lo_inner;
use crate::double::SemiDouble;
use crate::traits::Int as _;

pub(crate) fn atanh<F: Ln>(x: F) -> F {
    let e = x.raw_exp();
    let absx = x.abs();
    if (e == F::MAX_RAW_EXP && x.raw_mant() != F::Raw::ZERO)
        || e <= (F::EXP_OFFSET - F::RawExp::from(F::MANT_BITS))
    {
        // propagate NaN
        // or
        // very small, includes subnormal and zero
        // atanh(x) ~= x
        // also handles atanh(-0) = -0
        x
    } else if absx == F::ONE {
        // atanh(±1) = ±inf
        F::INFINITY.copysign(x)
    } else if x.abs() > F::ONE {
        // |x| > 1, return NaN
        F::NAN
    } else {
        atanh_inner(x)
    }
}

fn atanh_inner<F: Ln>(x: F) -> F {
    // t1 = 2 * x / (1 - x)
    let t1 = SemiDouble::new(F::TWO * x) / SemiDouble::new_qsub11(F::ONE, x);

    // t2 = (1 + x) / (1 - x) = t1 + 1
    let t2 = t1 + F::ONE;
    let t2 = t2.to_norm();

    // atanh(x) = 0.5 * ln((1 + x) / (1 - x))
    F::HALF * ln_hi_lo_inner(t2.hi(), t2.lo())
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::Float;

    fn test<F: Float + FloatMath>() {
        use crate::atanh;

        let f = F::parse;

        assert_is_nan!(atanh(F::NAN));
        assert_is_nan!(atanh(f("1.5")));
        assert_is_nan!(atanh(f("-1.5")));
        assert_is_nan!(atanh(F::INFINITY));
        assert_is_nan!(atanh(F::NEG_INFINITY));
        assert_total_eq!(atanh(F::ZERO), F::ZERO);
        assert_total_eq!(atanh(-F::ZERO), -F::ZERO);
        assert_total_eq!(atanh(F::ONE), F::INFINITY);
        assert_total_eq!(atanh(-F::ONE), F::NEG_INFINITY);
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
