use crate::traits::{Float, Int as _};

pub(crate) trait Trigonometric: Float {
    fn sin_finite(x: Self) -> Self;

    fn cos_finite(x: Self) -> Self;

    fn sin_cos_finite(x: Self) -> (Self, Self);

    fn tan_finite(x: Self) -> Self;

    fn sind_finite(x: Self) -> Self;

    fn cosd_finite(x: Self) -> Self;

    fn sind_cosd_finite(x: Self) -> (Self, Self);

    fn tand_finite(x: Self) -> Self;

    fn sinpi_finite(x: Self) -> Self;

    fn cospi_finite(x: Self) -> Self;

    fn sinpi_cospi_finite(x: Self) -> (Self, Self);

    fn tanpi_finite(x: Self) -> Self;
}

pub(crate) fn sin<F: Trigonometric>(x: F) -> F {
    let xexp = x.raw_exp();
    if xexp == F::MAX_RAW_EXP {
        // sin(NaN) = NaN
        // sin(±inf) = NaN
        F::NAN
    } else if xexp == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // sin(±0) = ±0
        x
    } else {
        F::sin_finite(x)
    }
}

pub(crate) fn cos<F: Trigonometric>(x: F) -> F {
    let xexp = x.raw_exp();
    if xexp == F::MAX_RAW_EXP {
        // cos(NaN) = NaN
        // cos(±inf) = NaN
        F::NAN
    } else if xexp == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // cos(±0) = ±1
        F::ONE
    } else {
        F::cos_finite(x)
    }
}

pub(crate) fn sin_cos<F: Trigonometric>(x: F) -> (F, F) {
    let xexp = x.raw_exp();
    if xexp == F::MAX_RAW_EXP {
        // sin(NaN) = NaN
        // sin(±inf) = NaN
        // cos(NaN) = NaN
        // cos(±inf) = NaN
        (F::NAN, F::NAN)
    } else if xexp == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // sin(±0) = ±0
        // cos(±0) = 1
        (x, F::ONE)
    } else {
        F::sin_cos_finite(x)
    }
}

pub(crate) fn tan<F: Trigonometric>(x: F) -> F {
    let xexp = x.raw_exp();
    if xexp == F::MAX_RAW_EXP {
        // tan(NaN) = NaN
        // tan(±inf) = NaN
        F::NAN
    } else if xexp == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // tan(±0) = ±0
        x
    } else {
        F::tan_finite(x)
    }
}

pub(crate) fn sind<F: Trigonometric>(x: F) -> F {
    let xexp = x.raw_exp();
    if xexp == F::MAX_RAW_EXP {
        // sind(NaN) = NaN
        // sind(±inf) = NaN
        F::NAN
    } else if xexp == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // sind(±0) = ±0
        x
    } else {
        F::sind_finite(x)
    }
}

pub(crate) fn cosd<F: Trigonometric>(x: F) -> F {
    let xexp = x.raw_exp();
    if xexp == F::MAX_RAW_EXP {
        // cosd(NaN) = NaN
        // cosd(±inf) = NaN
        F::NAN
    } else if xexp == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // cosd(±0) = ±1
        F::ONE
    } else {
        F::cosd_finite(x)
    }
}

pub(crate) fn sind_cosd<F: Trigonometric>(x: F) -> (F, F) {
    let xexp = x.raw_exp();
    if xexp == F::MAX_RAW_EXP {
        // sind(NaN) = NaN
        // sind(±inf) = NaN
        // cosd(NaN) = NaN
        // cosd(±inf) = NaN
        (F::NAN, F::NAN)
    } else if xexp == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // sind(±0) = ±0
        // cosd(±0) = 1
        (x, F::ONE)
    } else {
        F::sind_cosd_finite(x)
    }
}

pub(crate) fn tand<F: Trigonometric>(x: F) -> F {
    let xexp = x.raw_exp();
    if xexp == F::MAX_RAW_EXP {
        // tand(NaN) = NaN
        // tand(±inf) = NaN
        F::NAN
    } else if xexp == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // tand(±0) = ±0
        x
    } else {
        F::tand_finite(x)
    }
}

pub(crate) fn sinpi<F: Trigonometric>(x: F) -> F {
    let xexp = x.raw_exp();
    if xexp == F::MAX_RAW_EXP {
        // sinpi(NaN) = NaN
        // sinpi(±inf) = NaN
        F::NAN
    } else if xexp == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // sinpi(±0) = ±0
        x
    } else {
        F::sinpi_finite(x)
    }
}

pub(crate) fn cospi<F: Trigonometric>(x: F) -> F {
    let xexp = x.raw_exp();
    if xexp == F::MAX_RAW_EXP {
        // cospi(NaN) = NaN
        // cospi(±inf) = NaN
        F::NAN
    } else if xexp == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // cospi(±0) = ±1
        F::ONE
    } else {
        F::cospi_finite(x)
    }
}

pub(crate) fn sinpi_cospi<F: Trigonometric>(x: F) -> (F, F) {
    let xexp = x.raw_exp();
    if xexp == F::MAX_RAW_EXP {
        // sinpi(NaN) = NaN
        // sinpi(±inf) = NaN
        // cospi(NaN) = NaN
        // cospi(±inf) = NaN
        (F::NAN, F::NAN)
    } else if xexp == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // sinpi(±0) = ±0
        // cospi(±0) = 1
        (x, F::ONE)
    } else {
        F::sinpi_cospi_finite(x)
    }
}

pub(crate) fn tanpi<F: Trigonometric>(x: F) -> F {
    let xexp = x.raw_exp();
    if xexp == F::MAX_RAW_EXP {
        // tanpi(NaN) = NaN
        // tanpi(±inf) = NaN
        F::NAN
    } else if xexp == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // tanpi(±0) = ±0
        x
    } else {
        F::tanpi_finite(x)
    }
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::Float;

    fn test_sin_cos<F: Float + FloatMath>() {
        use crate::{cos, sin, sin_cos};

        let test_nan = |arg: F| {
            let sin1 = sin(arg);
            let cos1 = cos(arg);
            let (sin2, cos2) = sin_cos(arg);
            assert_is_nan!(sin1);
            assert_is_nan!(cos1);
            assert_is_nan!(sin2);
            assert_is_nan!(cos2);
        };

        let test_value = |arg: F, expected_sin: F, expected_cos: F| {
            let sin1 = sin(arg);
            let cos1 = cos(arg);
            let (sin2, cos2) = sin_cos(arg);
            assert_total_eq!(sin1, expected_sin);
            assert_total_eq!(cos1, expected_cos);
            assert_total_eq!(sin2, expected_sin);
            assert_total_eq!(cos2, expected_cos);
        };

        test_nan(F::NAN);
        test_nan(F::INFINITY);
        test_nan(F::NEG_INFINITY);
        test_value(F::ZERO, F::ZERO, F::ONE);
        test_value(-F::ZERO, -F::ZERO, F::ONE);
    }

    fn test_tan<F: Float + FloatMath>() {
        use crate::tan;

        assert_is_nan!(tan(F::NAN));
        assert_is_nan!(tan(F::INFINITY));
        assert_is_nan!(tan(F::NEG_INFINITY));
        assert_total_eq!(tan(F::ZERO), F::ZERO);
        assert_total_eq!(tan(-F::ZERO), -F::ZERO);
    }

    fn test_sind_cosd<F: Float + FloatMath>() {
        use crate::{cosd, sind, sind_cosd};

        let f = F::parse;

        let test_nan = |arg: F| {
            let sin1 = sind(arg);
            let cos1 = cosd(arg);
            let (sin2, cos2) = sind_cosd(arg);
            assert_is_nan!(sin1);
            assert_is_nan!(cos1);
            assert_is_nan!(sin2);
            assert_is_nan!(cos2);
        };

        let test_value = |arg: F, expected_sin: F, expected_cos: F| {
            let sin1 = sind(arg);
            let cos1 = cosd(arg);
            let (sin2, cos2) = sind_cosd(arg);
            assert_total_eq!(sin1, expected_sin);
            assert_total_eq!(cos1, expected_cos);
            assert_total_eq!(sin2, expected_sin);
            assert_total_eq!(cos2, expected_cos);
        };

        test_nan(F::NAN);
        test_nan(F::INFINITY);
        test_nan(F::NEG_INFINITY);
        test_value(F::ZERO, F::ZERO, F::ONE);
        test_value(-F::ZERO, -F::ZERO, F::ONE);
        test_value(f("90"), F::ONE, F::ZERO);
        test_value(f("-90"), -F::ONE, F::ZERO);
        test_value(f("180"), F::ZERO, -F::ONE);
        test_value(f("-180"), -F::ZERO, -F::ONE);
        test_value(f("270"), -F::ONE, F::ZERO);
        test_value(f("-270"), F::ONE, F::ZERO);
        test_value(f("360"), F::ZERO, F::ONE);
        test_value(f("-360"), -F::ZERO, F::ONE);
        test_value(f("450"), F::ONE, F::ZERO);
        test_value(f("-450"), -F::ONE, F::ZERO);
    }

    fn test_tand<F: Float + FloatMath>() {
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

    fn test_sinpi_cospi<F: Float + FloatMath>() {
        use crate::{cospi, sinpi, sinpi_cospi};

        let f = F::parse;

        let test_nan = |arg: F| {
            let sin1 = sinpi(arg);
            let cos1 = cospi(arg);
            let (sin2, cos2) = sinpi_cospi(arg);
            assert_is_nan!(sin1);
            assert_is_nan!(cos1);
            assert_is_nan!(sin2);
            assert_is_nan!(cos2);
        };

        let test_value = |arg: F, expected_sin: F, expected_cos: F| {
            let sin1 = sinpi(arg);
            let cos1 = cospi(arg);
            let (sin2, cos2) = sinpi_cospi(arg);
            assert_total_eq!(sin1, expected_sin);
            assert_total_eq!(cos1, expected_cos);
            assert_total_eq!(sin2, expected_sin);
            assert_total_eq!(cos2, expected_cos);
        };

        test_nan(F::NAN);
        test_nan(F::INFINITY);
        test_nan(F::NEG_INFINITY);
        test_value(F::ZERO, F::ZERO, F::ONE);
        test_value(-F::ZERO, -F::ZERO, F::ONE);
        test_value(f("0.5"), F::ONE, F::ZERO);
        test_value(f("-0.5"), -F::ONE, F::ZERO);
        test_value(f("1.0"), F::ZERO, -F::ONE);
        test_value(f("-1.0"), -F::ZERO, -F::ONE);
        test_value(f("1.5"), -F::ONE, F::ZERO);
        test_value(f("-1.5"), F::ONE, F::ZERO);
        test_value(f("2.0"), F::ZERO, F::ONE);
        test_value(f("-2.0"), -F::ZERO, F::ONE);
        test_value(f("2.5"), F::ONE, F::ZERO);
        test_value(f("-2.5"), -F::ONE, F::ZERO);
    }

    fn test_tanpi<F: Float + FloatMath>() {
        use crate::tanpi;

        let f = F::parse;

        assert_is_nan!(tanpi(F::NAN));
        assert_is_nan!(tanpi(F::INFINITY));
        assert_is_nan!(tanpi(F::NEG_INFINITY));
        assert_total_eq!(tanpi(F::ZERO), F::ZERO);
        assert_total_eq!(tanpi(-F::ZERO), -F::ZERO);
        assert_total_eq!(tanpi(f("0.5")), F::INFINITY);
        assert_total_eq!(tanpi(f("-0.5")), F::NEG_INFINITY);
        assert_total_eq!(tanpi(f("1.0")), -F::ZERO);
        assert_total_eq!(tanpi(f("-1.0")), F::ZERO);
        assert_total_eq!(tanpi(f("1.5")), F::NEG_INFINITY);
        assert_total_eq!(tanpi(f("-1.5")), F::INFINITY);
        assert_total_eq!(tanpi(f("2.0")), F::ZERO);
        assert_total_eq!(tanpi(f("-2.0")), -F::ZERO);
        assert_total_eq!(tanpi(f("2.5")), F::INFINITY);
        assert_total_eq!(tanpi(f("-2.5")), F::NEG_INFINITY);
    }

    #[test]
    fn test_f32() {
        test_sin_cos::<f32>();
        test_tan::<f32>();
        test_sind_cosd::<f32>();
        test_tand::<f32>();
        test_sinpi_cospi::<f32>();
        test_tanpi::<f32>();
    }

    #[test]
    fn test_f64() {
        test_sin_cos::<f64>();
        test_tan::<f64>();
        test_sind_cosd::<f64>();
        test_tand::<f64>();
        test_sinpi_cospi::<f64>();
        test_tanpi::<f64>();
    }
}
