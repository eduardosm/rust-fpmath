use super::{round_as_i_f, scalbn_medium};
use crate::double::DenormDouble;
use crate::traits::{CastInto as _, Float, Int as _, Like};

pub(crate) trait Exp<L = Like<Self>>: Float {
    fn log2_e() -> Self;
    fn ln_2_hi() -> Self;
    fn ln_2_lo() -> Self;
    fn exp_lo_th() -> Self;
    fn exp_hi_th() -> Self;
    fn exp_m1_lo_th() -> Self;
    fn exp_m1_hi_th() -> Self;

    fn exp_special_poly(x2: Self) -> Self;

    fn exp_m1_special_poly(x2: Self) -> Self;
}

/// Calculates `e^x`
pub(crate) fn exp<F: Exp>(x: F) -> F {
    if x >= F::exp_hi_th() {
        // also handles x = inf
        F::INFINITY
    } else if x <= F::exp_lo_th() {
        // also handles x = -inf
        F::ZERO
    } else {
        let e = x.raw_exp();
        if e == F::RawExp::ZERO {
            // x is zero or subnormal
            // exp(x) ~= 1
            F::one()
        } else if e == F::MAX_RAW_EXP {
            // x is NaN, propagate
            x
        } else {
            exp_inner(x)
        }
    }
}

pub(crate) fn exp_m1<F: Exp>(x: F) -> F {
    if x >= F::exp_m1_hi_th() {
        // also handles x = inf
        F::INFINITY
    } else if x <= F::exp_m1_lo_th() {
        // also handles x = -inf
        -F::one()
    } else {
        let e = x.raw_exp();
        if e == F::RawExp::ZERO || e == F::MAX_RAW_EXP {
            // x is zero or subnormal: exp(x) - 1 ~= x
            // or
            // x is NaN: propagate
            x
        } else {
            exp_m1_inner(x)
        }
    }
}

pub(super) fn exp_inner<F: Exp>(x: F) -> F {
    // Split x into k, r_hi, r_lo such as:
    //  - x = k*ln(2) + r_hi + r_lo
    //  - k is an integer
    //  - |r| <= 0.5*ln(2)
    let (k, r_hi, r_lo) = exp_split(x);

    // exp(x) = exp(k*ln(2) + r_hi + r_lo)
    exp_inner_common(k, r_hi, r_lo)
}

/// Calculates `exp(k*ln(2) + r_hi + r_lo)`
pub(super) fn exp_inner_common<F: Exp>(k: i32, r_hi: F, r_lo: F) -> F {
    // Based on the algorithm used by the msun math library

    let r = r_hi + r_lo;
    let r2 = r * r;

    // t1 = 2 - 2 * r / (exp(r) - 1)
    let t1 = r + F::exp_special_poly(r2);

    // t2 = exp(r) = 1 + r + (r * t1) / (2 - t1)
    //             = 1 + r_hi + r_lo + (r * t1) / (2 - t1)
    let t2 = F::one() + (r_hi + (r_lo + r * t1 / (F::two() - t1)));

    // exp(x) = exp(r_hi + r_lo) * 2^k = t2 * 2^k
    scalbn_medium(t2, k)
}

/// Calculates `exp(r_hi + r_lo)`
pub(super) fn hi_lo_exp_inner_common<F: Exp>(r_hi: F, r_lo: F) -> DenormDouble<F> {
    // Based on the algorithm used by the msun math library

    let r = r_hi + r_lo;
    let r2 = r * r;

    // t1 = 2 - 2 * r / (exp(r) - 1)
    let t1 = DenormDouble::new_qadd11(r, F::exp_special_poly(r2));

    // t2 = (r * t1) / (2 - t1)
    let rt1 = DenormDouble::new(r_hi, r_lo).to_semi() * t1.to_semi();
    let twomt1 = t1.qrsub1(F::two());
    let t2 = rt1.to_semi() / twomt1.to_semi();

    // t3 = exp(r) = 1 + r + t2
    let t3 = DenormDouble::new(r_hi, r_lo).qadd2(t2).qradd1(F::one());

    // exp(x) = exp(r_hi + r_lo) * 2^k = t2 * 2^k
    DenormDouble::new(t3.hi(), t3.lo())
}

fn exp_m1_inner<F: Exp>(x: F) -> F {
    // Based on the algorithm used by the msun math library

    // pseudo-consts
    let three = F::one() + F::two();
    let six = three * F::two();

    // Split x into k, r_hi, r_lo such as:
    //  - x = k*ln(2) + r_hi + r_lo
    //  - k is an integer
    //  - |r_hi| <= 0.5*ln(2)
    let (k, r_hi, r_lo) = exp_split(x);
    let (r_hi, r_lo) = F::norm_hi_lo_full(r_hi, r_lo);

    let r2 = r_hi * r_hi;
    let hr = F::half() * r_hi;
    let hr2 = F::half() * r2;

    // t1 = 6/r * ((exp(r) + 1) / (exp(r) - 1) - 2/r)
    let t1 = F::exp_m1_special_poly(r2);
    // t2 = 3 - t1 * 0.5 * r
    let t2 = three - t1 * hr;
    // t3 = 0.5 * r^2 * (t1 - t2) / (6 - r * t2)
    let t3 = hr2 * ((t1 - t2) / (six - r_hi * t2));
    // t4 = exp(r_hi + r_lo) - 1 - r_hi = r_hi * (r_lo - t3) + r_lo + 0.5 * r^2
    let t4 = (r_hi * (r_lo - t3) + r_lo) + hr2;

    // exp(x) - 1 = exp(r_hi + r_lo) * 2^k - 1 = (r_hi + t4 + 1) * 2^k - 1
    if k < F::MAX_EXP.into() {
        let s1 = F::exp2i_fast(k.cast_into());
        let sr = r_hi * s1;
        let st4 = t4 * s1;

        let t5 = DenormDouble::new_qadd11(s1, sr);
        let t6 = t5.qadd1(st4);
        let t7 = t6.qsub1(F::one());

        t7.to_single()
    } else {
        scalbn_medium((r_hi + t4) + F::one(), k)
    }
}

/// Splits `x` into `(k, r_hi, r_lo)`
///
/// Such as:
/// * `x = k*ln(2) + r_hi + r_lo`
/// * `k` is an integer
/// * `|r| <= 0.5*ln(2)`
pub(super) fn exp_split<F: Exp>(x: F) -> (i32, F, F) {
    let y = x * F::log2_e();
    let (k, kf) = round_as_i_f(y);
    // `kf * LN_2_HI` is exact because the lower bits of `LN_2_HI` are zero.
    let r_hi = x - kf * F::ln_2_hi();
    let r_lo = -kf * F::ln_2_lo();

    (k, r_hi, r_lo)
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test_exp<F: Float + FloatMath>(lo_th: &str, hi_th: &str) {
        use crate::exp;

        let lo_th = F::parse(lo_th);
        let hi_th = F::parse(hi_th);

        assert_is_nan!(exp(F::NAN));
        assert_total_eq!(exp(F::INFINITY), F::INFINITY);
        assert_total_eq!(exp(F::neg_infinity()), F::ZERO);
        assert_total_eq!(exp(F::ZERO), F::one());
        assert_total_eq!(exp(-F::ZERO), F::one());
        assert_total_eq!(exp(lo_th), F::ZERO);
        assert_total_eq!(exp(lo_th - F::one()), F::ZERO);
        assert_total_eq!(exp(lo_th - F::two()), F::ZERO);
        assert_total_eq!(exp(hi_th), F::INFINITY);
        assert_total_eq!(exp(hi_th + F::one()), F::INFINITY);
        assert_total_eq!(exp(hi_th + F::two()), F::INFINITY);
    }

    fn test_exp_m1<F: Float + FloatMath>(lo_th: &str, hi_th: &str) {
        use crate::exp_m1;

        let lo_th = F::parse(lo_th);
        let hi_th = F::parse(hi_th);

        assert_is_nan!(exp_m1(F::NAN));
        assert_total_eq!(exp_m1(F::INFINITY), F::INFINITY);
        assert_total_eq!(exp_m1(F::neg_infinity()), -F::one());
        assert_total_eq!(exp_m1(F::ZERO), F::ZERO);
        assert_total_eq!(exp_m1(-F::ZERO), -F::ZERO);
        assert_total_eq!(exp_m1(lo_th), -F::one());
        assert_total_eq!(exp_m1(lo_th - F::one()), -F::one());
        assert_total_eq!(exp_m1(lo_th - F::two()), -F::one());
        assert_total_eq!(exp_m1(hi_th), F::INFINITY);
        assert_total_eq!(exp_m1(hi_th + F::one()), F::INFINITY);
        assert_total_eq!(exp_m1(hi_th + F::two()), F::INFINITY);
    }

    #[test]
    fn test_f32() {
        test_exp::<f32>("-103.99", "88.9");
        test_exp_m1::<f32>("-87.9", "88.9");
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f32() {
        test_exp::<crate::SoftF32>("-103.99", "88.9");
        test_exp_m1::<crate::SoftF32>("-87.9", "88.9");
    }

    #[test]
    fn test_f64() {
        test_exp::<f64>("-745.9", "709.9");
        test_exp_m1::<f64>("-708.9", "709.9");
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f64() {
        test_exp::<crate::SoftF64>("-745.9", "709.9");
        test_exp_m1::<crate::SoftF64>("-708.9", "709.9");
    }
}
