use crate::traits::{Float, Int as _};

pub(crate) fn frexp<F: Float>(x: F) -> (F, i32) {
    let (y, edelta) = x.normalize_arg();
    let yexp = y.raw_exp();
    if yexp == F::RawExp::ZERO {
        // zero
        (F::ZERO.copysign(x), 0)
    } else if yexp == F::MAX_RAW_EXP {
        // infinity or NaN
        (y, 0)
    } else {
        // finite
        (
            y.set_exp(-F::Exp::ONE),
            (F::raw_exp_to_exp(yexp) + F::Exp::ONE + edelta).into(),
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test<F: Float + FloatMath>() {
        use crate::frexp;

        let f = F::parse;

        let test_nan = || {
            let (mant, exp) = frexp(F::NAN);
            assert_is_nan!(mant);
            assert_eq!(exp, 0);
        };

        let test = |x: F, expected_mant: F, expected_exp: i32| {
            let (mant, exp) = frexp(x);
            assert_total_eq!(mant, expected_mant);
            assert_eq!(exp, expected_exp);
        };

        test_nan();
        test(F::INFINITY, F::INFINITY, 0);
        test(F::neg_infinity(), F::neg_infinity(), 0);
        test(F::ZERO, F::ZERO, 0);
        test(-F::ZERO, -F::ZERO, 0);
        test(f("0.09375"), f("0.75"), -3);
        test(f("-0.09375"), f("-0.75"), -3);
        test(f("0.25"), f("0.5"), -1);
        test(f("-0.25"), f("-0.5"), -1);
        test(f("0.5"), f("0.5"), 0);
        test(f("-0.5"), f("-0.5"), 0);
        test(F::one(), f("0.5"), 1);
        test(-F::one(), f("-0.5"), 1);
        test(f("20"), f("0.625"), 5);
        test(f("-20"), f("-0.625"), 5);

        // Subnormal numbers
        let min_subnormal_exp: i32 =
            (F::MIN_NORMAL_EXP - F::Exp::try_from(F::MANT_BITS).ok().unwrap()).into();

        let r = |raw: u8| F::from_raw(raw.into());

        test(r(0b01), f("0.5"), min_subnormal_exp + 1);
        test(-r(0b01), f("-0.5"), min_subnormal_exp + 1);
        test(r(0b10), f("0.5"), min_subnormal_exp + 2);
        test(-r(0b10), f("-0.5"), min_subnormal_exp + 2);
        test(r(0b11), f("0.75"), min_subnormal_exp + 2);
        test(-r(0b11), f("-0.75"), min_subnormal_exp + 2);
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
