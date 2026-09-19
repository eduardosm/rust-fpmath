use super::{Reduce90Deg, Tan, reduce_90_deg, tan::tan_quadrant};

pub(crate) fn tand<F: Reduce90Deg + Tan>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        // tand(inf or NaN) = NaN
        F::NAN
    } else if e <= F::RawExp::from(F::MANT_BITS) {
        // very small, includes subnormal and zero
        // tand(x) ~= x * (π/180)
        // also handles tand(-0) = -0
        x * F::DEG_TO_RAD
    } else {
        let (n, y) = reduce_90_deg(x);
        tan_quadrant(n, x.sign(), y.hi(), y.lo())
    }
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::Float;

    fn test<F: Float + FloatMath>() {
        use crate::tand;

        let f = F::parse;

        assert_is_nan!(tand(F::NAN));
        assert_is_nan!(tand(F::INFINITY));
        assert_is_nan!(tand(F::NEG_INFINITY));
        assert_total_eq!(tand(F::ZERO), F::ZERO);
        assert_total_eq!(tand(-F::ZERO), -F::ZERO);
        assert_total_eq!(tand(f("90")), F::INFINITY);
        assert_total_eq!(tand(f("-90")), F::NEG_INFINITY);
        assert_total_eq!(tand(f("180")), -F::ZERO);
        assert_total_eq!(tand(f("-180")), F::ZERO);
        assert_total_eq!(tand(f("270")), F::NEG_INFINITY);
        assert_total_eq!(tand(f("-270")), F::INFINITY);
        assert_total_eq!(tand(f("360")), F::ZERO);
        assert_total_eq!(tand(f("-360")), -F::ZERO);
        assert_total_eq!(tand(f("450")), F::INFINITY);
        assert_total_eq!(tand(f("-450")), F::NEG_INFINITY);
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
