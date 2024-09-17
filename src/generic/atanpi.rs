use super::atan::{atan2_inner, atan_inner};
use super::{Atan, DivPi};
use crate::traits::{CastFrom as _, CastInto as _, Int as _};

pub(crate) fn atanpi<F: Atan + DivPi>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        if x.raw_mant() == F::Raw::ZERO {
            // atanpi(±inf) = ±0.5
            F::half().copysign(x)
        } else {
            // propagate NaN
            x
        }
    } else if e == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // atanpi(±0) = ±0
        x
    } else if e <= F::RawExp::from(F::MANT_BITS) {
        // very small, atanpi(x) ~= x / π

        // scale temporarily to avoid temporary subnormal numbers
        let logscale = F::Exp::TWO * F::Exp::cast_from(F::MANT_BITS);
        let scale = F::exp2i_fast(logscale);
        let descale = F::exp2i_fast(-logscale);

        let (x_hi, x_lo) = (x * scale).split_hi_lo();
        let y_hi = x_hi * F::frac_1_pi_hi();
        let y_lo = x_hi * F::frac_1_pi_lo() + x_lo * F::frac_1_pi();
        (y_hi + y_lo) * descale
    } else {
        let (y_hi, y_lo) = atan_inner(x);
        let (y_hi, y_lo) = F::norm_hi_lo_splitted(y_hi, y_lo);

        y_hi * F::frac_1_pi_hi() + (y_hi * F::frac_1_pi_lo() + y_lo * F::frac_1_pi())
    }
}

pub(crate) fn atan2pi<F: Atan + DivPi>(y: F, x: F) -> F {
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
        let quarter = F::half() * F::half();
        let three_quarter = quarter + quarter + quarter;
        // x = ±inf, y = ±inf
        match (nx.sign(), ny.sign()) {
            (false, false) => quarter,
            (false, true) => -quarter,
            (true, false) => three_quarter,
            (true, true) => -three_quarter,
        }
    } else if nxexp == F::MAX_RAW_EXP {
        // x = ±inf
        if nx.sign() {
            F::one().copysign(ny)
        } else {
            F::ZERO.copysign(ny)
        }
    } else if nyexp == F::MAX_RAW_EXP {
        // y = ±inf
        F::half().copysign(ny)
    } else if nyexp == F::RawExp::ZERO {
        // y = ±0
        if nx.sign() {
            F::one().copysign(ny)
        } else {
            ny
        }
    } else if nxexp == F::RawExp::ZERO {
        // x = ±0
        F::half().copysign(ny)
    } else if !nx.sign()
        && nxexp > nyexp
        && (nxexp - nyexp) >= ((F::MAX_RAW_EXP >> 1) - F::MANT_BITS.into())
    {
        let scale = F::exp2i_fast(F::Exp::cast_from(F::MANT_BITS));
        let descale = F::exp2i_fast(-F::Exp::cast_from(F::MANT_BITS));

        // y/x is very small
        // atan2pi(y, x) ~= (y/x) / π
        let (ny_hi, ny_lo) = (ny * scale).split_hi_lo();
        let (nx_hi, nx_lo) = nx.split_hi_lo();

        let (nydeg_hi, nydeg_lo) = (ny_hi * F::frac_1_pi_hi()).split_hi_lo();
        let nydeg_lo = nydeg_lo + (ny_hi * F::frac_1_pi_lo() + ny_lo * F::frac_1_pi());

        let (t1_hi, t1_lo) = F::div_hi_lo(nydeg_hi, nydeg_lo, nx_hi, nx_lo);
        (t1_hi + t1_lo) * descale
    } else {
        let (y_hi, y_lo) = atan2_inner(ny, nx);
        let (y_hi, y_lo) = F::norm_hi_lo_splitted(y_hi, y_lo);

        y_hi * F::frac_1_pi_hi() + (y_hi * F::frac_1_pi_lo() + y_lo * F::frac_1_pi())
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test_atanpi<F: Float + FloatMath>() {
        use crate::atanpi;

        assert_is_nan!(atanpi(F::NAN));
        assert_total_eq!(atanpi(F::INFINITY), F::half());
        assert_total_eq!(atanpi(F::neg_infinity()), -F::half());
        assert_total_eq!(atanpi(F::ZERO), F::ZERO);
        assert_total_eq!(atanpi(-F::ZERO), -F::ZERO);
    }

    fn test_atan2pi<F: Float + FloatMath>() {
        use crate::atan2pi;

        let f = F::parse;

        assert_is_nan!(atan2pi(F::NAN, F::one()));
        assert_is_nan!(atan2pi(F::NAN, F::ZERO));
        assert_is_nan!(atan2pi(F::NAN, F::INFINITY));
        assert_is_nan!(atan2pi(F::NAN, F::NAN));
        assert_is_nan!(atan2pi(F::INFINITY, F::NAN));
        assert_is_nan!(atan2pi(F::ZERO, F::NAN));
        assert_is_nan!(atan2pi(F::one(), F::NAN));
        assert_total_eq!(atan2pi(F::ZERO, F::ZERO), F::ZERO);
        assert_total_eq!(atan2pi(-F::ZERO, F::ZERO), -F::ZERO);
        assert_total_eq!(atan2pi(F::ZERO, F::one()), F::ZERO);
        assert_total_eq!(atan2pi(-F::ZERO, F::one()), -F::ZERO);
        assert_total_eq!(atan2pi(F::ZERO, F::INFINITY), F::ZERO);
        assert_total_eq!(atan2pi(-F::ZERO, F::INFINITY), -F::ZERO);
        assert_total_eq!(atan2pi(F::ZERO, -F::ZERO), F::one());
        assert_total_eq!(atan2pi(-F::ZERO, -F::ZERO), -F::one());
        assert_total_eq!(atan2pi(F::ZERO, -F::one()), F::one());
        assert_total_eq!(atan2pi(-F::ZERO, -F::one()), -F::one());
        assert_total_eq!(atan2pi(F::INFINITY, F::ZERO), F::half());
        assert_total_eq!(atan2pi(F::INFINITY, -F::ZERO), F::half());
        assert_total_eq!(atan2pi(F::INFINITY, F::one()), F::half());
        assert_total_eq!(atan2pi(F::INFINITY, -F::one()), F::half());
        assert_total_eq!(atan2pi(F::neg_infinity(), F::ZERO), -F::half());
        assert_total_eq!(atan2pi(F::neg_infinity(), -F::ZERO), -F::half());
        assert_total_eq!(atan2pi(F::neg_infinity(), F::one()), -F::half());
        assert_total_eq!(atan2pi(F::neg_infinity(), -F::one()), -F::half());
        assert_total_eq!(atan2pi(F::ZERO, F::INFINITY), F::ZERO);
        assert_total_eq!(atan2pi(-F::ZERO, F::INFINITY), -F::ZERO);
        assert_total_eq!(atan2pi(F::one(), F::INFINITY), F::ZERO);
        assert_total_eq!(atan2pi(-F::one(), F::INFINITY), -F::ZERO);
        assert_total_eq!(atan2pi(F::ZERO, F::neg_infinity()), F::one());
        assert_total_eq!(atan2pi(-F::ZERO, F::neg_infinity()), -F::one());
        assert_total_eq!(atan2pi(F::one(), F::neg_infinity()), F::one());
        assert_total_eq!(atan2pi(-F::one(), F::neg_infinity()), -F::one());
        assert_total_eq!(atan2pi(F::INFINITY, F::INFINITY), f("0.25"));
        assert_total_eq!(atan2pi(F::neg_infinity(), F::INFINITY), f("-0.25"));
        assert_total_eq!(atan2pi(F::INFINITY, F::neg_infinity()), f("0.75"));
        assert_total_eq!(atan2pi(F::neg_infinity(), F::neg_infinity()), f("-0.75"));
    }

    #[test]
    fn test_f32() {
        test_atanpi::<f32>();
        test_atan2pi::<f32>();
    }

    #[test]
    fn test_soft_f32() {
        test_atanpi::<crate::SoftF32>();
        test_atan2pi::<crate::SoftF32>();
    }

    #[test]
    fn test_f64() {
        test_atanpi::<f64>();
        test_atan2pi::<f64>();
    }

    #[test]
    fn test_soft_f64() {
        test_atanpi::<crate::SoftF64>();
        test_atan2pi::<crate::SoftF64>();
    }
}
