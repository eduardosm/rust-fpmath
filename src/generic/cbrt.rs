use crate::traits::{Float, Int as _, Like};

pub(crate) trait Cbrt<L = Like<Self>>: Float {
    fn cbrt_2_hi() -> Self;
    fn cbrt_2_lo() -> Self;
    fn cbrt_4_hi() -> Self;
    fn cbrt_4_lo() -> Self;

    fn exp_mod_3(e: Self::Exp) -> i8;

    fn inv_cbrt_poly(x: Self) -> Self;
}

/// Returns `(sign, k, r)` as needed by `cbrt_inner`
fn cbrt_split<F: Cbrt>(x: F, edelta: F::Exp) -> (bool, F::Exp, F, F, F) {
    let k = x.exponent() + edelta;
    // 0 <= kmod3 <= 2
    let kmod3 = F::exp_mod_3(k);

    let (cb0_hi, cb0_lo) = match kmod3 {
        0 => (F::one(), F::ZERO),
        1 => (F::cbrt_2_hi(), F::cbrt_2_lo()),
        2 => (F::cbrt_4_hi(), F::cbrt_4_lo()),
        _ => unreachable!(),
    };

    // 1 <= r < 2
    let r = x.abs().set_exp(F::Exp::ZERO);

    (x.sign(), k - F::Exp::from(kmod3), r, cb0_hi, cb0_lo)
}

fn cbrt_inner<F: Cbrt>(x: F, edelta: F::Exp) -> F {
    let inv_three = F::one() / (F::one() + F::two());

    // Split x * 2^edelta = (-1)^sign * 2^k * r * cb0^3 such as
    // * k is an integer
    // * k mod 3 = 0
    // * 1 <= r < 2
    // * cb0 = cb0_hi + cb0_lo = 1, cbrt(2) or cbrt(4)
    let (sign, k, r, cb0_hi, cb0_lo) = cbrt_split(x, edelta);

    // Based on the algorithm used in SLEEF.
    // https://github.com/shibatch/sleef/wiki/Divisionless-iterative-approximation-method-of-cube-root

    // ta ~= cbrt(1 / r) with a polynomial approximation
    let ta = F::inv_cbrt_poly(r);

    // tb ~= cbrt(1 / r) with a Newton iteration
    // tb = ta - (1 / 3) * (r * ta^4 - ta)
    let ta2 = ta * ta;
    let ta4 = ta2 * ta2;
    let tb = (ta - inv_three * (r * ta4 - ta)).purify();

    // cbrt(r) = r * cbrt(1 / r)^2
    // cbrt(1 / r) is calculated with another Newton iteration
    let (tb2_hi, tb2_lo) = (tb * tb).split_hi_lo();

    let (tb4_hi, tb4_lo) = (tb2_hi * tb2_hi).split_hi_lo();
    let tb4_lo = tb4_lo + (F::two() * tb2_hi * tb2_lo + tb2_lo * tb2_lo);

    // tb4r = tb^4 * r
    let (r_hi, r_lo) = r.split_hi_lo();
    let tb4r_hi = r_hi * tb4_hi;
    let tb4r_lo = r * tb4_lo + r_lo * tb4_hi;

    // tc = r * tb^4 - tb
    let tc = (tb4r_hi - tb) + tb4r_lo;

    // td = (-2 / 3) * tb * (r * tb^4 - tb) = (-2 / 3) * tb * tc
    let td = ((-F::two() * inv_three) * tb * tc).purify();

    // te = tb^2 + (-2 / 3) * tb * (r * tb^4 - tb) = tb^2 + td
    let te_hi = (tb2_hi + td).split_hi();
    let te_lo = tb2_lo + ((tb2_hi - te_hi) + td);

    // tf = te * r = cbrt(r)
    let (tf_hi, tf_lo) = (te_hi * r_hi).split_hi_lo();
    let tf_lo = tf_lo + (te_hi * r_lo + te_lo * r);

    // tg = cbrt(r) * cb0 = tf * cb0
    let tg = tf_hi * cb0_hi + (tf_hi * cb0_lo + tf_lo * (cb0_hi + cb0_lo));

    // y = cbrt(r) * 2^(k / 3) * cb0 = tg * 2^(k / 3)
    let y = tg * F::exp2i_fast(k / (F::Exp::ONE + F::Exp::TWO));
    // cbrt(x) = (-1)^sign * cbrt(r) * 2^(k / 3) * cb0 = (-1)^sign * y
    F::from_raw(y.to_raw() | (F::Raw::from(sign) << (F::BITS - 1)))
}

pub(crate) fn cbrt<F: Cbrt>(x: F) -> F {
    let (y, edelta) = x.normalize_arg();
    let yexp = y.raw_exp();
    if yexp == F::RawExp::ZERO || yexp == F::MAX_RAW_EXP {
        // cbrt(±0) = ±0
        // or
        // propagate infinity or NaN
        y
    } else {
        cbrt_inner(y, edelta)
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test<F: Float + FloatMath>(
        full_e_mants: impl Clone + Iterator<Item = u64>,
        extra_e: impl Iterator<Item = i32>,
        extra_e_mants: impl Clone + Iterator<Item = u64>,
    ) {
        use crate::{cbrt, scalbn};

        assert_is_nan!(cbrt(F::NAN));
        assert_total_eq!(cbrt(F::neg_infinity()), F::neg_infinity());
        assert_total_eq!(cbrt(F::INFINITY), F::INFINITY);
        assert_total_eq!(cbrt(F::ZERO), F::ZERO);
        assert_total_eq!(cbrt(-F::ZERO), -F::ZERO);

        let min_normal_exp: i32 = F::MIN_NORMAL_EXP.into();
        let e_limit = (-min_normal_exp) / 3;
        for e in (-e_limit)..e_limit {
            for m in full_e_mants.clone() {
                let x = scalbn(F::cast_from(m), e - m.ilog2() as i32);
                let x3 = x * x * x;
                assert_total_eq!(cbrt(x3), x);
                assert_total_eq!(cbrt(-x3), -x);
            }
        }
        for e in extra_e {
            for m in extra_e_mants.clone() {
                let x = scalbn(F::cast_from(m), e - m.ilog2() as i32);
                let x3 = x * x * x;
                assert_total_eq!(cbrt(x3), x);
                assert_total_eq!(cbrt(-x3), -x);
            }
        }
    }

    #[test]
    fn test_f32() {
        test::<f32>(0x80..=0xFF, core::iter::empty(), core::iter::empty());
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f32() {
        test::<crate::SoftF32>(0x80..=0xFF, core::iter::empty(), core::iter::empty());
    }

    #[test]
    fn test_f64() {
        test::<f64>(0x100..=0x1FF, [-340, 0, 340].into_iter(), 0x10000..=0x1FFFF);
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f64() {
        test::<crate::SoftF64>(0x100..=0x1FF, [-340, 0, 340].into_iter(), 0x10000..=0x1FFFF);
    }
}
