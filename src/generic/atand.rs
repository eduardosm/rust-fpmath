use super::atan::{atan2_inner, atan_inner};
use super::{Atan, RadToDeg};
use crate::double::SemiDouble;
use crate::traits::{CastFrom as _, CastInto as _, Int as _};

pub(crate) fn atand<F: Atan + RadToDeg>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        if x.raw_mant() == F::Raw::ZERO {
            // atand(±inf) = ±90
            F::cast_from(90u32).copysign(x)
        } else {
            // propagate NaN
            x
        }
    } else if e <= F::MANT_BITS.into() {
        // subnormal or zero
        // atand(x) ~= x * (180/π)
        // also handles atand(-0) = -0
        x * F::rad_to_deg()
    } else {
        let y = atan_inner(x).to_semi();

        (y * F::rad_to_deg_ex()).to_single()
    }
}

pub(crate) fn atan2d<F: Atan + RadToDeg>(y: F, x: F) -> F {
    let (ny, nx) = if y.raw_exp() <= F::MANT_BITS.into() || x.raw_exp() <= F::MANT_BITS.into() {
        // convert possible subnormals to normals
        let scale = F::exp2i_fast((F::MANT_BITS * 2 + 1).cast_into());
        (y * scale, x * scale)
    } else {
        (y, x)
    };

    let nxexp = nx.raw_exp();
    let nyexp = ny.raw_exp();
    if (nxexp == F::MAX_RAW_EXP && nx.raw_mant() != F::Raw::ZERO)
        || (nyexp == F::MAX_RAW_EXP && ny.raw_mant() != F::Raw::ZERO)
    {
        // x and/or y is NaN
        F::NAN
    } else if nxexp == F::MAX_RAW_EXP && nyexp == F::MAX_RAW_EXP {
        // x = ±inf, y = ±inf
        match (nx.sign(), ny.sign()) {
            (false, false) => F::cast_from(45u32),
            (false, true) => -F::cast_from(45u32),
            (true, false) => F::cast_from(135u32),
            (true, true) => -F::cast_from(135u32),
        }
    } else if nxexp == F::MAX_RAW_EXP {
        // x = ±inf
        if nx.sign() {
            F::cast_from(180u32).copysign(ny)
        } else {
            F::ZERO.copysign(ny)
        }
    } else if nyexp == F::MAX_RAW_EXP {
        // y = ±inf
        F::cast_from(90u32).copysign(ny)
    } else if nyexp == F::RawExp::ZERO {
        // y = ±0
        if nx.sign() {
            F::cast_from(180u32).copysign(ny)
        } else {
            ny
        }
    } else if nxexp == F::RawExp::ZERO {
        // x = ±0
        F::cast_from(90u32).copysign(ny)
    } else if !nx.sign()
        && nxexp > nyexp
        && (nxexp - nyexp) >= ((F::MAX_RAW_EXP >> 1) - F::MANT_BITS.into())
    {
        let scale = F::exp2i_fast(F::Exp::cast_from(F::MANT_BITS));
        let descale = F::exp2i_fast(-F::Exp::cast_from(F::MANT_BITS));

        // y/x is very small
        // atan2d(y, x) ~= (y/x) * (180/π)
        let ny = SemiDouble::new(ny * scale);
        let nx = SemiDouble::new(nx);

        let nydeg = ny * F::rad_to_deg_ex();

        (nydeg.to_semi() / nx).to_single() * descale
    } else {
        let y = atan2_inner(ny, nx).to_semi();

        (y * F::rad_to_deg_ex()).to_single()
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test_atand<F: Float + FloatMath>() {
        use crate::atand;

        let f = F::parse;

        assert_is_nan!(atand(F::NAN));
        assert_total_eq!(atand(F::INFINITY), f("90"));
        assert_total_eq!(atand(F::neg_infinity()), f("-90"));
        assert_total_eq!(atand(F::ZERO), F::ZERO);
        assert_total_eq!(atand(-F::ZERO), -F::ZERO);
    }

    fn test_atan2d<F: Float + FloatMath>() {
        use crate::atan2d;

        let f = F::parse;

        assert_is_nan!(atan2d(F::NAN, F::one()));
        assert_is_nan!(atan2d(F::NAN, F::ZERO));
        assert_is_nan!(atan2d(F::NAN, F::INFINITY));
        assert_is_nan!(atan2d(F::NAN, F::NAN));
        assert_is_nan!(atan2d(F::INFINITY, F::NAN));
        assert_is_nan!(atan2d(F::ZERO, F::NAN));
        assert_is_nan!(atan2d(F::one(), F::NAN));
        assert_total_eq!(atan2d(F::ZERO, F::ZERO), F::ZERO);
        assert_total_eq!(atan2d(-F::ZERO, F::ZERO), -F::ZERO);
        assert_total_eq!(atan2d(F::ZERO, F::one()), F::ZERO);
        assert_total_eq!(atan2d(-F::ZERO, F::one()), -F::ZERO);
        assert_total_eq!(atan2d(F::ZERO, F::INFINITY), F::ZERO);
        assert_total_eq!(atan2d(-F::ZERO, F::INFINITY), -F::ZERO);
        assert_total_eq!(atan2d(F::ZERO, -F::ZERO), f("180"));
        assert_total_eq!(atan2d(-F::ZERO, -F::ZERO), f("-180"));
        assert_total_eq!(atan2d(F::ZERO, -F::one()), f("180"));
        assert_total_eq!(atan2d(-F::ZERO, -F::one()), f("-180"));
        assert_total_eq!(atan2d(F::INFINITY, F::ZERO), f("90"));
        assert_total_eq!(atan2d(F::INFINITY, -F::ZERO), f("90"));
        assert_total_eq!(atan2d(F::INFINITY, F::one()), f("90"));
        assert_total_eq!(atan2d(F::INFINITY, -F::one()), f("90"));
        assert_total_eq!(atan2d(F::neg_infinity(), F::ZERO), f("-90"));
        assert_total_eq!(atan2d(F::neg_infinity(), -F::ZERO), f("-90"));
        assert_total_eq!(atan2d(F::neg_infinity(), F::one()), f("-90"));
        assert_total_eq!(atan2d(F::neg_infinity(), -F::one()), f("-90"));
        assert_total_eq!(atan2d(F::ZERO, F::INFINITY), F::ZERO);
        assert_total_eq!(atan2d(-F::ZERO, F::INFINITY), -F::ZERO);
        assert_total_eq!(atan2d(F::one(), F::INFINITY), F::ZERO);
        assert_total_eq!(atan2d(-F::one(), F::INFINITY), -F::ZERO);
        assert_total_eq!(atan2d(F::ZERO, F::neg_infinity()), f("180"));
        assert_total_eq!(atan2d(-F::ZERO, F::neg_infinity()), f("-180"));
        assert_total_eq!(atan2d(F::one(), F::neg_infinity()), f("180"));
        assert_total_eq!(atan2d(-F::one(), F::neg_infinity()), f("-180"));
        assert_total_eq!(atan2d(F::INFINITY, F::INFINITY), f("45"));
        assert_total_eq!(atan2d(F::neg_infinity(), F::INFINITY), f("-45"));
        assert_total_eq!(atan2d(F::INFINITY, F::neg_infinity()), f("135"));
        assert_total_eq!(atan2d(F::neg_infinity(), F::neg_infinity()), f("-135"));
    }

    #[test]
    fn test_f32() {
        test_atand::<f32>();
        test_atan2d::<f32>();
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f32() {
        test_atand::<crate::SoftF32>();
        test_atan2d::<crate::SoftF32>();
    }

    #[test]
    fn test_f64() {
        test_atand::<f64>();
        test_atan2d::<f64>();
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f64() {
        test_atand::<crate::SoftF64>();
        test_atan2d::<crate::SoftF64>();
    }
}
