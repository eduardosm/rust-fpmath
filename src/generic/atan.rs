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
        let (y_hi, y_lo) = atan_inner(x);
        y_hi + y_lo
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
        let (y_hi, y_lo) = atan2_inner(ny, nx);
        y_hi + y_lo
    }
}

pub(super) fn atan_inner<F: Atan>(x: F) -> (F, F) {
    if x.abs() <= F::one() {
        let (x_hi, x_lo) = x.split_hi_lo();
        atan_inner_common(x_hi, x_lo)
    } else {
        let (x_hi, x_lo) = x.split_hi_lo();
        let (y_hi, y_lo) = F::recip_hi_lo(x_hi, x_lo);
        let (y_hi, y_lo) = F::norm_hi_lo_splitted(y_hi, y_lo);

        // t1 = atan(1 / x)
        let (t1_hi, t1_lo) = atan_inner_common(y_hi, y_lo);

        // t2 = atan(x) = ±pi/2 - atan(1 / x)
        let off_hi = F::frac_pi_2_hi().copysign(x);
        let off_lo = F::frac_pi_2_lo().copysign(x);

        let t2_hi = (off_hi - t1_hi).purify();
        let t2_lo = (off_lo - t1_lo) + ((off_hi - t2_hi) - t1_hi);

        (t2_hi, t2_lo)
    }
}

pub(super) fn atan2_inner<F: Atan>(mut n: F, mut d: F) -> (F, F) {
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

    let (n_hi, n_lo) = n.split_hi_lo();
    let (d_hi, d_lo) = d.split_hi_lo();

    // z = n/d
    let (z_hi, z_lo) = F::div_hi_lo(n_hi, n_lo, d_hi, d_lo);
    let (z_hi, z_lo) = F::norm_hi_lo_splitted(z_hi, z_lo);

    // t1 = atan(n/d)
    let (t1_hi, t1_lo) = atan_inner_common(z_hi, z_lo);

    // t2 = off * π/2
    let t2_hi = F::frac_pi_2_hi() * off;
    let t2_lo = F::frac_pi_2_lo() * off;

    // t3 = atan2(y, x) = atan(n/d) + off * π/2 = t1 + t2
    let t3_hi = (t1_hi + t2_hi).purify();
    let t3_lo = (t1_lo + t2_lo) + ((t2_hi - t3_hi) + t1_hi);

    (t3_hi, t3_lo)
}

pub(super) fn atan_inner_common<F: Atan>(x_hi: F, x_lo: F) -> (F, F) {
    let (x2_hi, x2_lo) = (x_hi * x_hi).split_hi_lo();
    let x2_lo = x2_lo + (F::two() * x_hi * x_lo + x_lo * x_lo);
    let x2 = x2_hi + x2_lo;

    // t1 = (atan(x) - x) / x^3 - k3
    let (k3, t1) = F::atan_poly(x2);

    // t2 = (atan(x) - x) / x^3
    let t2_hi = (t1 + k3).split_hi();
    let t2_lo = (k3 - t2_hi) + t1;

    let (x3_hi, x3_lo) = (x_hi * x2_hi).split_hi_lo();
    let x3_lo = x3_lo + (x_hi * x2_lo + x_lo * x2_hi + x_lo * x2_lo);

    // t3 = atan(x) - x
    let t3_hi = t2_hi * x3_hi;
    let t3_lo = t2_hi * x3_lo + t2_lo * x3_hi + t2_lo * x3_lo;

    // t4 = atan(x)
    let t4_hi = (t3_hi + x_hi).purify();
    let t4_lo = (x_lo + t3_lo) + ((x_hi - t4_hi) + t3_hi);

    (t4_hi, t4_lo)
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
