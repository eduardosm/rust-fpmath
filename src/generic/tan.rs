use super::{reduce_pi_2, ReducePi2};
use crate::traits::{Float, Int as _, Like};

pub(crate) trait Tan<L = Like<Self>>: Float {
    /// Calculates `tan(x) - x`
    ///
    /// Where:
    /// * `x2 = x^2`
    /// * `x3 = x^3`
    fn tan_poly(x2: Self, x3: Self) -> Self;
}

pub(super) fn tan_inner<F: Tan>(x_hi: F, x_lo: F, inv: bool) -> F {
    // let y = 0.5 * x
    // tan(x) = 2 * tan(y) / (1 - tan(y)^2)

    let y_hi = F::half() * x_hi;
    let y_lo = F::half() * x_lo;
    let y2 = y_hi * y_hi;
    let y3 = y_hi * y2;

    // ta = tan(y) - y
    let ta = F::tan_poly(y2, y3);

    // tb = tan(y) = ta + y
    let tb_hi = (ta + y_hi).split_hi();
    let tb_lo = y_lo + ((y_hi - tb_hi) + ta);

    let tb2_hi = tb_hi * tb_hi;
    let tb2_lo = F::two() * tb_hi * tb_lo + tb_lo * tb_lo;

    // tc = 2 * tan(y) = 2 * tb
    let tc_hi = F::two() * tb_hi;
    let tc_lo = F::two() * tb_lo;

    // td = 1 - tan(y)^2 = 1 - tb2
    let td_hi = (F::one() - tb2_hi).split_hi();
    let td_lo = ((F::one() - td_hi) - tb2_hi) - tb2_lo;

    // tan(x) = 2 * tan(y) / (1 - tan(y)^2) = tc / td
    // Calculate tan(x) or -1/tan(x) (tc / td or td / tc)

    let (n_hi, n_lo, d_hi, d_lo) = if inv {
        (-td_hi, -td_lo, tc_hi, tc_lo)
    } else {
        (tc_hi, tc_lo, td_hi, td_lo)
    };

    if d_hi == F::ZERO {
        F::INFINITY.copysign(n_hi)
    } else {
        let (q_hi, q_lo) = F::div_hi_lo(n_hi, n_lo, d_hi, d_lo);
        q_hi + q_lo
    }
}

pub(crate) fn tan<F: ReducePi2 + Tan>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        // tan(inf or NaN) = NaN
        F::NAN
    } else if e <= F::RawExp::ONE {
        // very small, includes subnormal and zero
        // tan(x) ~= x
        // also handles tan(-0) = -0
        x
    } else {
        let (n, y_hi, y_lo) = reduce_pi_2(x);

        tan_inner(y_hi, y_lo, (n & 1) != 0)
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test<F: Float + FloatMath>() {
        use crate::tan;

        assert_is_nan!(tan(F::NAN));
        assert_is_nan!(tan(F::INFINITY));
        assert_is_nan!(tan(F::neg_infinity()));
        assert_total_eq!(tan(F::ZERO), F::ZERO);
        assert_total_eq!(tan(-F::ZERO), -F::ZERO);
    }

    #[test]
    fn test_f32() {
        test::<f32>();
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f32() {
        test::<crate::SoftF32>();
    }

    #[test]
    fn test_f64() {
        test::<f64>();
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f64() {
        test::<crate::SoftF64>();
    }
}
