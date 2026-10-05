use crate::traits::{Float, Int as _};

pub(crate) trait Cbrt: Float {
    /// Calculates `cbrt(x)` for a finite and non-zero (but possibly
    /// subnormal) `x`.
    fn cbrt_finite(x: Self) -> Self;
}

pub(crate) fn cbrt<F: Cbrt>(x: F) -> F {
    let abs_raw = x.to_raw() & !F::SIGN_MASK;
    if abs_raw != F::Raw::ZERO && abs_raw < F::EXP_MASK {
        // x is finite and non-zero
        F::cbrt_finite(x)
    } else {
        // cbrt(±0) = ±0
        // or
        // propagate infinity or NaN
        x
    }
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::Float;

    fn test<F: Float + FloatMath>() {
        use crate::cbrt;

        assert_is_nan!(cbrt(F::NAN));
        assert_total_eq!(cbrt(F::NEG_INFINITY), F::NEG_INFINITY);
        assert_total_eq!(cbrt(F::INFINITY), F::INFINITY);
        assert_total_eq!(cbrt(F::ZERO), F::ZERO);
        assert_total_eq!(cbrt(-F::ZERO), -F::ZERO);
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
