use crate::traits::{CastInto as _, Float, Int as _};

pub(crate) fn round<F: Float>(x: F) -> F {
    let e = x.raw_exp();
    if e < (F::EXP_OFFSET - F::RawExp::ONE) {
        // abs(x) < 0.5
        // return zero without losing the sign
        F::ZERO.copysign(x)
    } else if e < F::EXP_OFFSET {
        // 0.5 <= abs(x) < 1
        // return ±1 keeping the sign
        F::one().copysign(x)
    } else {
        // x is NaN or abs(x) >= 1 (including infinity)
        // split integer and fractional parts
        // when NaN, infinity or exp >= MANT_BITS, fmask = 0
        let fmask = F::MANT_MASK >> (e - F::EXP_OFFSET).min(F::RawExp::from(F::MANT_BITS));
        let xraw = x.to_raw();
        let fpart = xraw & fmask;
        let ipart = xraw & !fmask;
        // add 1 to integer part if frac >= 0.5
        if fpart > (fmask / F::Raw::TWO) {
            F::from_raw(ipart + fmask + F::Raw::ONE)
        } else {
            F::from_raw(ipart)
        }
    }
}

/// Returns `x` rounded to the nearest integer as both integer and float.
///
/// `x` must be finite and `abs(int) < 2^min(31, MANT_BITS)`
pub(crate) fn round_as_i_f<F: Float>(x: F) -> (i32, F) {
    let e = x.raw_exp();
    if e < (F::EXP_OFFSET - F::RawExp::ONE) {
        // abs(x) < 0.5
        (0, F::ZERO)
    } else if e < F::EXP_OFFSET {
        // 0.5 <= abs(x) < 1
        (1 - (i32::from(x.sign()) << 1), F::one().copysign(x))
    } else {
        // 1 <= abs(x) < 2^min(31, MANT_BITS)
        let shift = F::RawExp::from(F::MANT_BITS) - (e - F::EXP_OFFSET);
        let imask = F::Raw::MAX << shift;
        let fmask = !imask;
        let xraw = x.to_raw();
        let fpart = xraw & fmask;
        let mut ipart_raw = xraw & !fmask;
        let mut ipart_i: i32 = (x.mant() >> shift).cast_into();
        if fpart > (fmask / F::Raw::TWO) {
            // frac >= 0.5
            ipart_raw += fmask + F::Raw::ONE;
            ipart_i += 1;
        }
        let ipart_f = F::from_raw(ipart_raw);
        if x.sign() {
            ipart_i = -ipart_i;
        }
        (ipart_i, ipart_f)
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test_round<F: Float + FloatMath>() {
        use crate::round;

        let one = F::one();
        let pt_1 = F::parse("0.1");
        let pt_5 = F::parse("0.5");
        let pt_9 = F::parse("0.9");

        assert_is_nan!(round(F::NAN));
        assert_total_eq!(round(F::INFINITY), F::INFINITY);
        assert_total_eq!(round(F::neg_infinity()), F::neg_infinity());

        for i in 0..20u32 {
            let x = F::cast_from(i);

            assert_total_eq!(round(x), x);
            assert_total_eq!(round(-x), -x);
            assert_total_eq!(round(x + pt_1), x);
            assert_total_eq!(round(-(x + pt_1)), -x);
            assert_total_eq!(round(x + pt_5), x + one);
            assert_total_eq!(round(-(x + pt_5)), -(x + one));
            assert_total_eq!(round(x + pt_9), x + one);
            assert_total_eq!(round(-(x + pt_9)), -(x + one));
        }
    }

    fn test_round_as_i_f<F: Float>() {
        let test = |x: F| {
            let (ipart_i, ipart_f) = super::round_as_i_f(x);
            let fpart = x - ipart_f;
            assert!(fpart.abs() <= F::half());
            assert_eq!(ipart_f, F::cast_from(ipart_i));
            assert_eq!(fpart + ipart_f, x);
        };

        let one_eight = F::parse("0.125");

        for i in 0..=1000u32 {
            for f in 0..8u32 {
                let x = F::cast_from(i) + F::cast_from(f) * one_eight;
                test(x);
                test(-x);
            }
        }
    }

    #[test]
    fn test_f32() {
        test_round::<f32>();
        test_round_as_i_f::<f32>();
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f32() {
        test_round::<crate::SoftF32>();
        test_round_as_i_f::<crate::SoftF32>();
    }

    #[test]
    fn test_f64() {
        test_round::<f64>();
        test_round_as_i_f::<f64>();
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f64() {
        test_round::<crate::SoftF64>();
        test_round_as_i_f::<crate::SoftF64>();
    }
}
