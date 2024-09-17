use crate::traits::Float;

pub(crate) fn trunc<F: Float>(x: F) -> F {
    let e = x.raw_exp();
    if e < F::EXP_OFFSET {
        // abs(x) < 1
        // return zero without losing the sign
        F::ZERO.copysign(x)
    } else {
        // x is NaN or abs(x) >= 1 (including infinity)
        // mask away the fractional digits from the mantissa.
        // fmask = 0 when NaN, infinity or exp >= MANT_BITS
        let fmask = F::MANT_MASK >> (e - F::EXP_OFFSET).min(F::RawExp::from(F::MANT_BITS));
        F::from_raw(x.to_raw() & !fmask)
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test<F: Float + FloatMath>() {
        use crate::trunc;

        let pt_1 = F::parse("0.1");
        let pt_5 = F::parse("0.5");
        let pt_9 = F::parse("0.9");

        assert_is_nan!(trunc(F::NAN));
        assert_total_eq!(trunc(F::INFINITY), F::INFINITY);
        assert_total_eq!(trunc(F::neg_infinity()), F::neg_infinity());

        for i in 0..20u32 {
            let x = F::cast_from(i);

            assert_total_eq!(trunc(x), x);
            assert_total_eq!(trunc(-x), -x);
            assert_total_eq!(trunc(x + pt_1), x);
            assert_total_eq!(trunc(-(x + pt_1)), -x);
            assert_total_eq!(trunc(x + pt_5), x);
            assert_total_eq!(trunc(-(x + pt_5)), -x);
            assert_total_eq!(trunc(x + pt_9), x);
            assert_total_eq!(trunc(-(x + pt_9)), -x);
        }
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
