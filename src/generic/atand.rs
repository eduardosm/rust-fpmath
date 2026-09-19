use super::atan::{atan_inner, atan2_inner};
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
        x * F::RAD_TO_DEG
    } else {
        let y = atan_inner(x).to_semi();

        (y * F::RAD_TO_DEG_EX).to_single()
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

        let nydeg = ny * F::RAD_TO_DEG_EX;

        // Copysign is to ensure the sign is preserved when the result is zero
        ((nydeg.to_semi() / nx).to_single() * descale).copysign(y)
    } else {
        let y = atan2_inner(ny, nx).to_semi();

        (y * F::RAD_TO_DEG_EX).to_single()
    }
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::Float;

    fn test_atand<F: Float + FloatMath>() {
        use crate::atand;

        let f = F::parse;

        assert_is_nan!(atand(F::NAN));
        assert_total_eq!(atand(F::INFINITY), f("90"));
        assert_total_eq!(atand(F::NEG_INFINITY), f("-90"));
        assert_total_eq!(atand(F::ZERO), F::ZERO);
        assert_total_eq!(atand(-F::ZERO), -F::ZERO);
    }

    fn test_atan2d<F: Float + FloatMath>() {
        use crate::{atan2d, scalbn};

        let f = F::parse;

        assert_is_nan!(atan2d(F::NAN, F::ONE));
        assert_is_nan!(atan2d(F::NAN, F::ZERO));
        assert_is_nan!(atan2d(F::NAN, F::INFINITY));
        assert_is_nan!(atan2d(F::NAN, F::NAN));
        assert_is_nan!(atan2d(F::INFINITY, F::NAN));
        assert_is_nan!(atan2d(F::ZERO, F::NAN));
        assert_is_nan!(atan2d(F::ONE, F::NAN));
        assert_total_eq!(atan2d(F::ZERO, F::ZERO), F::ZERO);
        assert_total_eq!(atan2d(-F::ZERO, F::ZERO), -F::ZERO);
        assert_total_eq!(atan2d(F::ZERO, F::ONE), F::ZERO);
        assert_total_eq!(atan2d(-F::ZERO, F::ONE), -F::ZERO);
        assert_total_eq!(atan2d(F::ZERO, F::INFINITY), F::ZERO);
        assert_total_eq!(atan2d(-F::ZERO, F::INFINITY), -F::ZERO);
        assert_total_eq!(atan2d(F::ZERO, -F::ZERO), f("180"));
        assert_total_eq!(atan2d(-F::ZERO, -F::ZERO), f("-180"));
        assert_total_eq!(atan2d(F::ZERO, -F::ONE), f("180"));
        assert_total_eq!(atan2d(-F::ZERO, -F::ONE), f("-180"));
        assert_total_eq!(atan2d(F::INFINITY, F::ZERO), f("90"));
        assert_total_eq!(atan2d(F::INFINITY, -F::ZERO), f("90"));
        assert_total_eq!(atan2d(F::INFINITY, F::ONE), f("90"));
        assert_total_eq!(atan2d(F::INFINITY, -F::ONE), f("90"));
        assert_total_eq!(atan2d(F::NEG_INFINITY, F::ZERO), f("-90"));
        assert_total_eq!(atan2d(F::NEG_INFINITY, -F::ZERO), f("-90"));
        assert_total_eq!(atan2d(F::NEG_INFINITY, F::ONE), f("-90"));
        assert_total_eq!(atan2d(F::NEG_INFINITY, -F::ONE), f("-90"));
        assert_total_eq!(atan2d(F::ZERO, F::INFINITY), F::ZERO);
        assert_total_eq!(atan2d(-F::ZERO, F::INFINITY), -F::ZERO);
        assert_total_eq!(atan2d(F::ONE, F::INFINITY), F::ZERO);
        assert_total_eq!(atan2d(-F::ONE, F::INFINITY), -F::ZERO);
        assert_total_eq!(atan2d(F::ZERO, F::NEG_INFINITY), f("180"));
        assert_total_eq!(atan2d(-F::ZERO, F::NEG_INFINITY), f("-180"));
        assert_total_eq!(atan2d(F::ONE, F::NEG_INFINITY), f("180"));
        assert_total_eq!(atan2d(-F::ONE, F::NEG_INFINITY), f("-180"));
        assert_total_eq!(atan2d(F::INFINITY, F::INFINITY), f("45"));
        assert_total_eq!(atan2d(F::NEG_INFINITY, F::INFINITY), f("-45"));
        assert_total_eq!(atan2d(F::INFINITY, F::NEG_INFINITY), f("135"));
        assert_total_eq!(atan2d(F::NEG_INFINITY, F::NEG_INFINITY), f("-135"));

        let small = scalbn(F::ONE, F::MIN_NORMAL_EXP.into() / 2);
        assert_total_eq!(atan2d(small, F::LARGEST), F::ZERO);
        assert_total_eq!(atan2d(-small, F::LARGEST), -F::ZERO);
        assert_total_eq!(atan2d(small, -F::LARGEST), f("180"));
        assert_total_eq!(atan2d(-small, -F::LARGEST), f("-180"));
    }

    #[test]
    fn test_f32() {
        test_atand::<f32>();
        test_atan2d::<f32>();
    }

    #[test]
    fn test_f64() {
        test_atand::<f64>();
        test_atan2d::<f64>();
    }
}
