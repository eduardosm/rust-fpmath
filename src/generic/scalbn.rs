use crate::traits::{CastInto as _, Float};

pub(crate) fn scalbn<F: Float>(x: F, y: i32) -> F {
    let (y1, y2, y3) = scalbn_split3::<F>(y);

    let p1 = F::exp2i_fast(y1);
    let p2 = F::exp2i_fast(y2);
    let p3 = F::exp2i_fast(y3);
    ((x * p1) * p2) * p3
}

#[inline]
pub(crate) fn scalbn_medium<F: Float>(x: F, y: i32) -> F {
    let (y1, y2) = scalbn_split2::<F>(y);

    let p1 = F::exp2i_fast(y1);
    let p2 = F::exp2i_fast(y2);
    (x * p1) * p2
}

#[inline]
fn scalbn_split2<F: Float>(e: i32) -> (F::Exp, F::Exp) {
    if e >= 0 {
        let e1 = e.min(F::MAX_EXP.into());
        let e2 = (e - e1).min(F::MAX_EXP.into());
        (e1.cast_into(), e2.cast_into())
    } else {
        // largest exponent (in abs) last to avoid subnormal double rounding
        let e2 = e.max(F::MIN_NORMAL_EXP.into());
        let e1 = (e - e2).max(F::MIN_NORMAL_EXP.into());
        (e1.cast_into(), e2.cast_into())
    }
}

#[inline]
fn scalbn_split3<F: Float>(e: i32) -> (F::Exp, F::Exp, F::Exp) {
    if e >= 0 {
        let e1 = e.min(F::MAX_EXP.into());
        let e2 = (e - e1).min(F::MAX_EXP.into());
        let e3 = ((e - e1) - e2).min(F::MAX_EXP.into());
        (e1.cast_into(), e2.cast_into(), e3.cast_into())
    } else {
        // largest exponent (in abs) last to avoid subnormal double rounding
        let e3 = e.max(F::MIN_NORMAL_EXP.into());
        let e2 = (e - e3).max(F::MIN_NORMAL_EXP.into());
        let e1 = ((e - e3) - e2).max(F::MIN_NORMAL_EXP.into());
        (e1.cast_into(), e2.cast_into(), e3.cast_into())
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::{Float, Int as _};
    use crate::FloatMath;

    fn test<F: Float + FloatMath>() {
        use crate::scalbn;

        let f = F::parse;

        assert_is_nan!(scalbn(F::NAN, 0));
        assert_is_nan!(scalbn(F::NAN, i32::MAX));
        assert_is_nan!(scalbn(F::NAN, i32::MIN));

        assert_total_eq!(scalbn(F::INFINITY, 0), F::INFINITY);
        assert_total_eq!(scalbn(F::INFINITY, i32::MAX), F::INFINITY);
        assert_total_eq!(scalbn(F::INFINITY, i32::MIN), F::INFINITY);

        assert_total_eq!(scalbn(F::neg_infinity(), 0), F::neg_infinity());
        assert_total_eq!(scalbn(F::neg_infinity(), i32::MAX), F::neg_infinity());
        assert_total_eq!(scalbn(F::neg_infinity(), i32::MIN), F::neg_infinity());

        assert_total_eq!(scalbn(F::ZERO, 0), F::ZERO);
        assert_total_eq!(scalbn(-F::ZERO, 0), -F::ZERO);
        assert_total_eq!(scalbn(F::ZERO, i32::MAX), F::ZERO);
        assert_total_eq!(scalbn(-F::ZERO, i32::MAX), -F::ZERO);
        assert_total_eq!(scalbn(F::ZERO, i32::MIN), F::ZERO);
        assert_total_eq!(scalbn(-F::ZERO, i32::MIN), -F::ZERO);

        assert_total_eq!(scalbn(F::one(), i32::MIN), F::ZERO);
        assert_total_eq!(scalbn(-F::one(), i32::MIN), -F::ZERO);
        assert_total_eq!(scalbn(F::one(), i32::MAX), F::INFINITY);
        assert_total_eq!(scalbn(-F::one(), i32::MAX), F::neg_infinity());
        assert_total_eq!(scalbn(f("10"), 2), f("40"));
        assert_total_eq!(scalbn(f("-10"), 2), f("-40"));
        assert_total_eq!(scalbn(f("2.5"), 3), f("20"));
        assert_total_eq!(scalbn(f("-2.5"), 3), f("-20"));

        let min_subnormal_exp = F::MIN_NORMAL_EXP - F::Exp::try_from(F::MANT_BITS).ok().unwrap();
        let from_subnormal_exp = |e: F::Exp| F::from_raw(F::Raw::ONE << (e - min_subnormal_exp));

        let min_subnormal = from_subnormal_exp(min_subnormal_exp);
        let max_subnormal = from_subnormal_exp(F::MIN_NORMAL_EXP - F::Exp::ONE);
        let min_normal = F::exp2i_fast(F::MIN_NORMAL_EXP);
        let max_normal = F::exp2i_fast(F::MAX_EXP);

        let max_norm_to_min_sub: i32 = (min_subnormal_exp - F::MAX_EXP).into();
        let min_norm_to_max_norm: i32 = (F::MAX_EXP - F::MIN_NORMAL_EXP).into();

        assert_total_eq!(scalbn(min_normal, 1), min_normal * F::two());
        assert_total_eq!(scalbn(-min_normal, 1), -min_normal * F::two());
        assert_total_eq!(scalbn(max_normal, max_norm_to_min_sub), min_subnormal);
        assert_total_eq!(scalbn(-max_normal, max_norm_to_min_sub), -min_subnormal);
        assert_total_eq!(scalbn(max_normal, max_norm_to_min_sub - 1), F::ZERO);
        assert_total_eq!(scalbn(-max_normal, max_norm_to_min_sub - 1), -F::ZERO);
        assert_total_eq!(scalbn(max_normal, 1), F::INFINITY);
        assert_total_eq!(scalbn(-max_normal, 1), F::neg_infinity());
        assert_total_eq!(scalbn(min_normal, min_norm_to_max_norm), max_normal);
        assert_total_eq!(scalbn(-min_normal, min_norm_to_max_norm), -max_normal);
        assert_total_eq!(scalbn(min_normal, min_norm_to_max_norm + 1), F::INFINITY);
        assert_total_eq!(
            scalbn(-min_normal, min_norm_to_max_norm + 1),
            F::neg_infinity()
        );

        let min_sub_to_max_exp: i32 = (F::MAX_EXP - min_subnormal_exp).into();

        let mant_bits_m1 = F::MANT_BITS - 1;
        assert_total_eq!(scalbn(min_subnormal, 1), min_subnormal * F::two());
        assert_total_eq!(scalbn(-min_subnormal, 1), -min_subnormal * F::two());
        assert_total_eq!(scalbn(min_subnormal, mant_bits_m1.into()), max_subnormal);
        assert_total_eq!(scalbn(-min_subnormal, mant_bits_m1.into()), -max_subnormal);
        assert_total_eq!(scalbn(min_subnormal, -1), F::ZERO);
        assert_total_eq!(scalbn(-min_subnormal, -1), -F::ZERO);
        assert_total_eq!(scalbn(min_subnormal, min_sub_to_max_exp + 1), F::INFINITY);
        assert_total_eq!(
            scalbn(-min_subnormal, min_sub_to_max_exp + 1),
            F::neg_infinity()
        );
        assert_total_eq!(scalbn(min_subnormal, min_sub_to_max_exp), max_normal);
        assert_total_eq!(scalbn(-min_subnormal, min_sub_to_max_exp), -max_normal);
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
