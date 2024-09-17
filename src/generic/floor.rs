use crate::traits::{Float, Int as _};

pub(crate) fn floor<F: Float>(x: F) -> F {
    let e = x.raw_exp();
    if e < F::EXP_OFFSET {
        // abs(x) < 1
        if !x.sign() || (x.to_raw() & (F::EXP_MASK | F::MANT_MASK)) == F::Raw::ZERO {
            // 0 <= x < 1
            // return zero without losing the sign
            F::ZERO.copysign(x)
        } else {
            // -1 < x < 0
            // return -1.0
            -F::one()
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
        if x.sign() && fpart != F::Raw::ZERO {
            F::from_raw(ipart + fmask + F::Raw::ONE)
        } else {
            F::from_raw(ipart)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test<F: Float + FloatMath>() {
        use crate::floor;

        let one = F::one();
        let pt_1 = F::parse("0.1");
        let pt_5 = F::parse("0.5");
        let pt_9 = F::parse("0.9");

        assert_is_nan!(floor(F::NAN));
        assert_total_eq!(floor(F::INFINITY), F::INFINITY);
        assert_total_eq!(floor(F::neg_infinity()), F::neg_infinity());

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
