use super::round_fi;
use crate::traits::{CastFrom as _, CastInto as _, Float, Int as _};

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

/// Reduces the angle argument `x` (in degrees), returning `(n, y)`
/// such as:
/// * `|y| <= 45`
/// * `0 <= n <= 3`
/// * `x = 360*M + 90*n + y`
/// * `M` is an integer
pub(crate) fn reduce_90_deg<F: Float>(x: F) -> (u8, F) {
    let xabs = x.abs();
    let xexp = x.exponent();
    if xabs <= F::cast_from(45u32) {
        (0, x)
    } else if xexp <= F::Exp::cast_from(F::MANT_BITS - 4) {
        let (f_n, n) = round_fi(x * (F::ONE / F::cast_from(90u32)));
        let y = x - f_n * F::cast_from(90u32);
        let n: u8 = n.cast_into();

        (n & 3, y)
    } else if xexp < F::Exp::cast_from(F::BITS - 1) {
        let xraw = x.to_raw();

        // ixint = int part of |x|
        // xfrac = frac part of |x|
        let (ixint, xfrac) = if xexp <= F::Exp::cast_from(F::MANT_BITS) {
            let shift: u8 = (F::Exp::cast_from(F::MANT_BITS) - xexp).cast_into();

            let fxint = F::from_raw(xraw & !F::SIGN_MASK & (F::Raw::MAX << shift));
            let ixint = x.mant() >> shift;
            let xfrac = xabs - fxint;

            (ixint, xfrac)
        } else {
            let shift: u8 = (xexp - F::Exp::cast_from(F::MANT_BITS)).cast_into();
            let ixint = x.mant() << shift;

            (ixint, F::ZERO)
        };

        // n = round(trunc(|x|) / 90)
        let n = (ixint + F::Raw::from(45u8)) / F::Raw::from(90u8);
        // rem = trunc(|x|) - n * 90
        let irem = F::SRaw::cast_from(ixint) - F::SRaw::cast_from(F::Raw::from(90u8) * n);
        let frem: F = irem.cast_into();

        let y = xfrac + frem;

        let n: u8 = n.cast_into();
        if x.is_sign_negative() {
            (n.wrapping_neg() & 3, -y)
        } else {
            (n & 3, y)
        }
    } else {
        // |x| = xm * 2^xe
        let xm = x.mant();
        let xe = xexp - F::Exp::cast_from(F::MANT_BITS);
        let xe: u16 = xe.cast_into();

        // EXP2_MOD45[i] = mod(2^(i + 3), 45)
        const EXP2_MOD45: [u8; 12] = [1, 2, 4, 8, 16, 32, 19, 38, 31, 17, 34, 23];

        // t = xm * mod(2^xe, 45)
        debug_assert!(xe > 3);
        let t = xm * F::Raw::from(EXP2_MOD45[usize::from((xe - 3) % 12)]);
        // rem45 = mod(|x|, 45)
        let rem45: u16 = (t % F::Raw::from(45u8)).cast_into();
        // rem360 = mod(|x|, 360) = mod(|x|, 45) * 8
        let rem360 = rem45 * 8;

        let mut n = (rem360 / 90) as u8;
        let mut rem90 = (rem360 % 90) as i16;
        if rem90 >= 45 {
            n += 1;
            rem90 -= 90;
        }

        let y = F::cast_from(rem90);

        if x.is_sign_negative() {
            (n.wrapping_neg() & 3, -y)
        } else {
            (n & 3, y)
        }
    }
}

/// Reduces the angle argument `x` (in radians/π), returning `(n, y)`
/// such as:
/// * `|y| <= 0.25`
/// * `0 <= n <= 3`
/// * `x = 2*M + 0.5*n + y`
/// * `M` is an integer
pub(crate) fn reduce_half_revs<F: Float>(x: F) -> (u8, F) {
    let xexp = x.exponent();
    if xexp < -F::Exp::TWO {
        // |x| < 0.25
        (0, x)
    } else if xexp == -F::Exp::TWO {
        // 0.25 <= abs(x) < 0.5
        let fpart = x - F::HALF.copysign(x);
        let n = if x.is_sign_negative() { 3 } else { 1 };

        (n, fpart)
    } else if xexp < F::Exp::cast_from(F::MANT_BITS) {
        // Split x = fpart + ipart * 0.5
        let shift: u8 = (F::Exp::cast_from(F::MANT_BITS) - (xexp + F::Exp::ONE)).cast_into();
        let fmask = !(F::Raw::MAX << shift);
        let xraw = x.to_raw();
        let fpart = xraw & fmask;
        let ipart = xraw & !fmask;
        let mut ipart_i = x.mant() >> shift;
        let mut fpart_f = x - F::from_raw(ipart);
        // Round to nearest
        if fpart > (fmask / F::Raw::TWO) {
            fpart_f = fpart_f - F::HALF.copysign(x);
            ipart_i += F::Raw::ONE;
        }

        let mut n: u8 = ipart_i.cast_into();
        if x.is_sign_negative() {
            n = n.wrapping_neg();
        }

        (n & 3, fpart_f)
    } else if xexp == F::Exp::cast_from(F::MANT_BITS) {
        // The lowest bit of the integer part is zero
        let n: u8 = ((x.to_raw() & F::Raw::ONE) << 1).cast_into();
        if x.is_sign_negative() {
            (n, F::ZERO)
        } else {
            (n.wrapping_neg() & 3, F::ZERO)
        }
    } else {
        // The two lowest bits of the integer part are zero
        (0, F::ZERO)
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
