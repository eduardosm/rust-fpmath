use crate::traits::{CastFrom as _, Float, Int as _};

pub(crate) fn round<F: Float>(x: F) -> F {
    let e = x.raw_exp();
    if e < (F::EXP_OFFSET - F::RawExp::ONE) {
        // abs(x) < 0.5
        // return zero without losing the sign
        F::ZERO.copysign(x)
    } else if e < F::EXP_OFFSET {
        // 0.5 <= abs(x) < 1
        // return ±1 keeping the sign
        F::ONE.copysign(x)
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

/// Returns `x` rounded to the nearest integer as both float and integer.
///
/// `x` must be finite and `abs(x) < 2^MANT_BITS`
pub(crate) fn round_fi<F: Float>(x: F) -> (F, F::SRaw) {
    let e = x.raw_exp();
    if e < (F::EXP_OFFSET - F::RawExp::ONE) {
        // abs(x) < 0.5
        (F::ZERO, F::SRaw::ZERO)
    } else if e < F::EXP_OFFSET {
        // 0.5 <= abs(x) < 1
        (
            F::ONE.copysign(x),
            F::SRaw::ONE - (F::SRaw::from(x.is_sign_negative()) << 1),
        )
    } else {
        // 1 <= abs(x) < 2^MANT_BITS
        let shift = F::RawExp::from(F::MANT_BITS) - (e - F::EXP_OFFSET);
        let imask = F::Raw::MAX << shift;
        let fmask = !imask;
        let xraw = x.to_raw();
        let fpart = xraw & fmask;
        let mut ipart_raw = xraw & !fmask;
        let mut ipart_i = F::SRaw::cast_from(x.mant() >> shift);
        if fpart > (fmask / F::Raw::TWO) {
            // frac >= 0.5
            ipart_raw += fmask + F::Raw::ONE;
            ipart_i += F::SRaw::ONE;
        }
        let ipart_f = F::from_raw(ipart_raw);
        if x.is_sign_negative() {
            ipart_i = -ipart_i;
        }
        (ipart_f, ipart_i)
    }
}

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

pub(crate) fn floor<F: Float>(x: F) -> F {
    let e = x.raw_exp();
    if e < F::EXP_OFFSET {
        // abs(x) < 1
        if !x.is_sign_negative() || (x.to_raw() & (F::EXP_MASK | F::MANT_MASK)) == F::Raw::ZERO {
            // 0 <= x < 1
            // return zero without losing the sign
            F::ZERO.copysign(x)
        } else {
            // -1 < x < 0
            // return -1.0
            -F::ONE
        }
    } else {
        // x is NaN or abs(x) >= 1 (including infinity)
        // split integer and fractional parts
        // when NaN, infinity or exp >= MANT_BITS, fmask = 0
        let fmask = F::MANT_MASK >> (e - F::EXP_OFFSET).min(F::RawExp::from(F::MANT_BITS));
        let xraw = x.to_raw();
        let fpart = xraw & fmask;
        let ipart = xraw & !fmask;
        // add 1 to integer part if x is negative and there
        // are non-zero fractional digits
        if x.is_sign_negative() && fpart != F::Raw::ZERO {
            F::from_raw(ipart + fmask + F::Raw::ONE)
        } else {
            F::from_raw(ipart)
        }
    }
}

pub(crate) fn ceil<F: Float>(x: F) -> F {
    let e = x.raw_exp();
    if e < F::EXP_OFFSET {
        // abs(x) < 1
        if !x.is_sign_negative() && (x.to_raw() & (F::EXP_MASK | F::MANT_MASK)) != F::Raw::ZERO {
            // 0 < x < 1
            // return 1.0
            F::ONE
        } else {
            // -1 < x <= 0
            // return zero without losing the sign
            F::ZERO.copysign(x)
        }
    } else {
        // x is NaN or abs(x) >= 1 (including infinity)
        // split integer and fractional parts
        // when NaN, infinity or exp >= MANT_BITS, fmask = 0
        let fmask = F::MANT_MASK >> (e - F::EXP_OFFSET).min(F::RawExp::from(F::MANT_BITS));
        let xraw = x.to_raw();
        let fpart = xraw & fmask;
        let ipart = xraw & !fmask;
        // add 1 to integer part if x is positive and there
        // are non-zero fractional digits
        if !x.is_sign_negative() && fpart != F::Raw::ZERO {
            F::from_raw(ipart + fmask + F::Raw::ONE)
        } else {
            F::from_raw(ipart)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::{CastInto as _, Float};

    fn test_round<F: Float + FloatMath>() {
        use crate::round;

        let one = F::ONE;
        let pt_1 = F::parse("0.1");
        let pt_5 = F::parse("0.5");
        let pt_9 = F::parse("0.9");

        assert_is_nan!(round(F::NAN));
        assert_total_eq!(round(F::INFINITY), F::INFINITY);
        assert_total_eq!(round(F::NEG_INFINITY), F::NEG_INFINITY);

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

    fn test_round_fi<F: Float>() {
        let test = |x: F| {
            let (ipart_f, ipart_i) = super::round_fi(x);
            let fpart = x - ipart_f;
            assert!(fpart.abs() <= F::HALF);
            assert_eq!(ipart_f, ipart_i.cast_into());
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

    fn test_trunc<F: Float + FloatMath>() {
        use crate::trunc;

        let pt_1 = F::parse("0.1");
        let pt_5 = F::parse("0.5");
        let pt_9 = F::parse("0.9");

        assert_is_nan!(trunc(F::NAN));
        assert_total_eq!(trunc(F::INFINITY), F::INFINITY);
        assert_total_eq!(trunc(F::NEG_INFINITY), F::NEG_INFINITY);

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

    fn test_floor<F: Float + FloatMath>() {
        use crate::floor;

        let one = F::ONE;
        let pt_1 = F::parse("0.1");
        let pt_5 = F::parse("0.5");
        let pt_9 = F::parse("0.9");

        assert_is_nan!(floor(F::NAN));
        assert_total_eq!(floor(F::INFINITY), F::INFINITY);
        assert_total_eq!(floor(F::NEG_INFINITY), F::NEG_INFINITY);

        for i in 0..20u32 {
            let x = F::cast_from(i);

            assert_total_eq!(floor(x), x);
            assert_total_eq!(floor(-x), -x);
            assert_total_eq!(floor(x + pt_1), x);
            assert_total_eq!(floor(-(x + pt_1)), -(x + one));
            assert_total_eq!(floor(x + pt_5), x);
            assert_total_eq!(floor(-(x + pt_5)), -(x + one));
            assert_total_eq!(floor(x + pt_9), x);
            assert_total_eq!(floor(-(x + pt_9)), -(x + one));
        }
    }

    fn test_ceil<F: Float + FloatMath>() {
        use crate::ceil;

        let one = F::ONE;
        let pt_1 = F::parse("0.1");
        let pt_5 = F::parse("0.5");
        let pt_9 = F::parse("0.9");

        assert_is_nan!(ceil(F::NAN));
        assert_total_eq!(ceil(F::INFINITY), F::INFINITY);
        assert_total_eq!(ceil(F::NEG_INFINITY), F::NEG_INFINITY);

        for i in 0..20u32 {
            let x = F::cast_from(i);

            assert_total_eq!(ceil(x), x);
            assert_total_eq!(ceil(-x), -x);
            assert_total_eq!(ceil(x + pt_1), x + one);
            assert_total_eq!(ceil(-(x + pt_1)), -x);
            assert_total_eq!(ceil(x + pt_5), x + one);
            assert_total_eq!(ceil(-(x + pt_5)), -x);
            assert_total_eq!(ceil(x + pt_9), x + one);
            assert_total_eq!(ceil(-(x + pt_9)), -x);
        }
    }

    #[test]
    fn test_f32() {
        test_round::<f32>();
        test_round_fi::<f32>();
        test_trunc::<f32>();
        test_floor::<f32>();
        test_ceil::<f32>();
    }

    #[test]
    fn test_f64() {
        test_round::<f64>();
        test_round_fi::<f64>();
        test_trunc::<f64>();
        test_floor::<f64>();
        test_ceil::<f64>();
    }
}
