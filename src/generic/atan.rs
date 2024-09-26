use crate::double::{DenormDouble, SemiDouble};
use crate::traits::{CastInto as _, FloatConsts, Int as _, Like};

pub(crate) trait Atan<L = Like<Self>>: FloatConsts {
    fn frac_pi_2_hi() -> Self;
    fn frac_pi_2_lo() -> Self;
    fn frac_3pi_4() -> Self;

    // Returns `(k3, (atan(x) - x) / x^3 - k3)`
    fn atan_poly(x2: Self) -> (Self, Self);
}

pub(crate) fn atan<F: Atan>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        if x.raw_mant() == F::Raw::ZERO {
            // atan(±inf) = ±π/2
            F::frac_pi_2().copysign(x)
        } else {
            // propagate NaN
            x
        }
    } else if e <= F::MANT_BITS.into() {
        // subnormal or zero
        // atan(x) ~= x
        // also handles atan(-0) = -0
        x
    } else {
        atan_inner(x).to_single()
    }
}

pub(crate) fn atan2<F: Atan>(y: F, x: F) -> F {
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
            (false, false) => F::frac_pi_4(),
            (false, true) => -F::frac_pi_4(),
            (true, false) => F::frac_3pi_4(),
            (true, true) => -F::frac_3pi_4(),
        }
    } else if nxexp == F::MAX_RAW_EXP {
        // x = ±inf
        if nx.sign() {
            F::pi().copysign(ny)
        } else {
            F::ZERO.copysign(ny)
        }
    } else if nyexp == F::MAX_RAW_EXP {
        // y = ±inf
        F::frac_pi_2().copysign(ny)
    } else if nyexp == F::RawExp::ZERO {
        // y = ±0
        if nx.sign() {
            F::pi().copysign(ny)
        } else {
            ny
        }
    } else if nxexp == F::RawExp::ZERO {
        // x = ±0
        F::frac_pi_2().copysign(ny)
    } else if !nx.sign()
        && nxexp > nyexp
        && (nxexp - nyexp) >= ((F::MAX_RAW_EXP >> 1) - F::MANT_BITS.into())
    {
        // y/x is very small
        // atan2(y, x) ~= y/x
        ny / nx
    } else {
        atan2_inner(ny, nx).to_single()
    }
}

pub(super) fn atan_inner<F: Atan>(x: F) -> DenormDouble<F> {
    if x.abs() <= F::one() {
        atan_inner_common(SemiDouble::new(x))
    } else {
        let y = DenormDouble::new_recip(x);

        // t1 = atan(1 / x)
        let t1 = atan_inner_common(y.to_semi());

        // atan(x) = ±pi/2 - atan(1 / x)
        let off = DenormDouble::new(F::frac_pi_2_hi().copysign(x), F::frac_pi_2_lo().copysign(x));
        off.qsub2(t1)
    }
}

pub(super) fn atan2_inner<F: Atan>(mut n: F, mut d: F) -> DenormDouble<F> {
    let ysgn = n.sign();
    let xsgn = d.sign();

    let mut off = F::ZERO;
    if xsgn {
        off = F::two().set_sign(ysgn);
    }
    if n.abs() > d.abs() {
        core::mem::swap(&mut n, &mut d);
        n = -n;
        off = off + F::one().set_sign(ysgn ^ xsgn);
    }

    // z = n/d
    let z = DenormDouble::new_div11(n, d);

    // t1 = atan(n/d)
    let t1 = atan_inner_common(z.to_semi());

    // t2 = off * π/2
    let t2 = DenormDouble::new(F::frac_pi_2_hi() * off, F::frac_pi_2_lo() * off);

    // atan2(y, x) = atan(n/d) + off * π/2 = t1 + t2
    t2.qadd2(t1)
}

pub(super) fn atan_inner_common<F: Atan>(x: SemiDouble<F>) -> DenormDouble<F> {
    let x2 = x.square();

    // t1 = (atan(x) - x) / x^3 - k3
    let (k3, t1) = F::atan_poly(x2.to_single());

    // t2 = (atan(x) - x) / x^3 = t1 + k3
    let t2 = SemiDouble::new_qadd11(k3, t1);

    let x3 = x * x2.to_semi();

    // t3 = atan(x) - x = t2 * x^3
    let t3 = x3.to_semi() * t2;

    // atan(x) = t3 + x
    x.to_denorm().qadd2(t3)
}

#[cfg(test)]
mod tests {
    use super::Atan;
    use crate::FloatMath;

    fn test_atan<F: Atan + FloatMath>() {
        use crate::atan;

        assert_is_nan!(atan(F::NAN));
        assert_total_eq!(atan(F::INFINITY), F::frac_pi_2());
        assert_total_eq!(atan(F::neg_infinity()), -F::frac_pi_2());
        assert_total_eq!(atan(F::ZERO), F::ZERO);
        assert_total_eq!(atan(-F::ZERO), -F::ZERO);
    }

    fn test_atan2<F: Atan + FloatMath>() {
        use crate::atan2;

        assert_is_nan!(atan2(F::NAN, F::one()));
        assert_is_nan!(atan2(F::NAN, F::ZERO));
        assert_is_nan!(atan2(F::NAN, F::INFINITY));
        assert_is_nan!(atan2(F::NAN, F::NAN));
        assert_is_nan!(atan2(F::INFINITY, F::NAN));
        assert_is_nan!(atan2(F::ZERO, F::NAN));
        assert_is_nan!(atan2(F::one(), F::NAN));
        assert_total_eq!(atan2(F::ZERO, F::ZERO), F::ZERO);
        assert_total_eq!(atan2(-F::ZERO, F::ZERO), -F::ZERO);
        assert_total_eq!(atan2(F::ZERO, F::one()), F::ZERO);
        assert_total_eq!(atan2(-F::ZERO, F::one()), -F::ZERO);
        assert_total_eq!(atan2(F::ZERO, F::INFINITY), F::ZERO);
        assert_total_eq!(atan2(-F::ZERO, F::INFINITY), -F::ZERO);
        assert_total_eq!(atan2(F::ZERO, -F::ZERO), F::pi());
        assert_total_eq!(atan2(-F::ZERO, -F::ZERO), -F::pi());
        assert_total_eq!(atan2(F::ZERO, -F::one()), F::pi());
        assert_total_eq!(atan2(-F::ZERO, -F::one()), -F::pi());
        assert_total_eq!(atan2(F::INFINITY, F::ZERO), F::frac_pi_2());
        assert_total_eq!(atan2(F::INFINITY, -F::ZERO), F::frac_pi_2());
        assert_total_eq!(atan2(F::INFINITY, F::one()), F::frac_pi_2());
        assert_total_eq!(atan2(F::INFINITY, -F::one()), F::frac_pi_2());
        assert_total_eq!(atan2(F::neg_infinity(), F::ZERO), -F::frac_pi_2());
        assert_total_eq!(atan2(F::neg_infinity(), -F::ZERO), -F::frac_pi_2());
        assert_total_eq!(atan2(F::neg_infinity(), F::one()), -F::frac_pi_2());
        assert_total_eq!(atan2(F::neg_infinity(), -F::one()), -F::frac_pi_2());
        assert_total_eq!(atan2(F::ZERO, F::INFINITY), F::ZERO);
        assert_total_eq!(atan2(-F::ZERO, F::INFINITY), -F::ZERO);
        assert_total_eq!(atan2(F::one(), F::INFINITY), F::ZERO);
        assert_total_eq!(atan2(-F::one(), F::INFINITY), -F::ZERO);
        assert_total_eq!(atan2(F::ZERO, F::neg_infinity()), F::pi());
        assert_total_eq!(atan2(-F::ZERO, F::neg_infinity()), -F::pi());
        assert_total_eq!(atan2(F::one(), F::neg_infinity()), F::pi());
        assert_total_eq!(atan2(-F::one(), F::neg_infinity()), -F::pi());
        assert_total_eq!(atan2(F::INFINITY, F::INFINITY), F::frac_pi_4());
        assert_total_eq!(atan2(F::neg_infinity(), F::INFINITY), -F::frac_pi_4());
        assert_total_eq!(atan2(F::INFINITY, F::neg_infinity()), F::frac_3pi_4());
        assert_total_eq!(
            atan2(F::neg_infinity(), F::neg_infinity()),
            -F::frac_3pi_4()
        );
    }

    #[test]
    fn test_f32() {
        test_atan::<f32>();
        test_atan2::<f32>();
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f32() {
        test_atan::<crate::SoftF32>();
        test_atan2::<crate::SoftF32>();
    }

    #[test]
    fn test_f64() {
        test_atan::<f64>();
        test_atan2::<f64>();
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f64() {
        test_atan::<crate::SoftF64>();
        test_atan2::<crate::SoftF64>();
    }
}
