use crate::traits::{Float, Int as _};

pub(crate) trait Cbrt: Float {
    fn cbrt_finite(x: Self, edelta: Self::Exp) -> Self;
}

pub(crate) fn cbrt<F: Cbrt>(x: F) -> F {
    let (y, edelta) = x.normalize_arg();
    let yexp = y.raw_exp();
    if yexp == F::RawExp::ZERO || yexp == F::MAX_RAW_EXP {
        // cbrt(±0) = ±0
        // or
        // propagate infinity or NaN
        y
    } else {
        F::cbrt_finite(y, edelta)
    }
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::Float;

    fn test<F: Float + FloatMath>(
        full_e_mants: impl Clone + Iterator<Item = u64>,
        extra_e: impl Iterator<Item = i32>,
        extra_e_mants: impl Clone + Iterator<Item = u64>,
    ) {
        use crate::{cbrt, scalbn};

        assert_is_nan!(cbrt(F::NAN));
        assert_total_eq!(cbrt(F::NEG_INFINITY), F::NEG_INFINITY);
        assert_total_eq!(cbrt(F::INFINITY), F::INFINITY);
        assert_total_eq!(cbrt(F::ZERO), F::ZERO);
        assert_total_eq!(cbrt(-F::ZERO), -F::ZERO);

        let min_normal_exp: i32 = F::MIN_NORMAL_EXP.into();
        let e_limit = (-min_normal_exp) / 3;
        for e in (-e_limit)..e_limit {
            for m in full_e_mants.clone() {
                let x = scalbn(F::cast_from(m), e - m.ilog2() as i32);
                let x3 = x * x * x;
                assert_total_eq!(cbrt(x3), x);
                assert_total_eq!(cbrt(-x3), -x);
            }
        }
        for e in extra_e {
            for m in extra_e_mants.clone() {
                let x = scalbn(F::cast_from(m), e - m.ilog2() as i32);
                let x3 = x * x * x;
                assert_total_eq!(cbrt(x3), x);
                assert_total_eq!(cbrt(-x3), -x);
            }
        }
    }

    #[test]
    fn test_f32() {
        test::<f32>(0x80..=0xFF, core::iter::empty(), core::iter::empty());
    }

    #[test]
    fn test_f64() {
        test::<f64>(0x100..=0x1FF, [-340, 0, 340].into_iter(), 0x10000..=0x1FFFF);
    }
}
