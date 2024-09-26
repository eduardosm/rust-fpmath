use super::{reduce_pi_2, ReducePi2};
use crate::double::{DenormDouble, SemiDouble};
use crate::traits::{Float, Int as _, Like};

pub(crate) trait Tan<L = Like<Self>>: Float {
    /// Calculates `tan(x) - x`
    ///
    /// Where:
    /// * `x2 = x^2`
    /// * `x3 = x^3`
    fn tan_poly(x2: Self, x3: Self) -> Self;
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

pub(super) fn tan_inner<F: Tan>(x_hi: F, x_lo: F, inv: bool) -> F {
    // let y = 0.5 * x
    // tan(x) = 2 * tan(y) / (1 - tan(y)^2)

    let y = DenormDouble::new(x_hi, x_lo).pmul1(F::half());
    let y2 = y.hi() * y.hi();
    let y3 = y.hi() * y2;

    // ta = tan(y) - y
    let ta = F::tan_poly(y2, y3);

    // tb = tan(y) = ta + y
    let tb = SemiDouble::new_qadd21(y, ta);
    let tb2 = tb.square();

    // tc = 2 * tan(y) = 2 * tb
    let tc = tb.pmul1(F::two());

    // td = 1 - tan(y)^2 = 1 - tb2
    let td = SemiDouble::new_qsub12(F::one(), tb2);

    // tan(x) = 2 * tan(y) / (1 - tan(y)^2) = tc / td
    // Calculate tan(x) or -1/tan(x) (tc / td or td / tc)

    let (n, d) = if inv { (-td, tc) } else { (tc, td) };

    if d.hi() == F::ZERO {
        F::INFINITY.copysign(n.hi())
    } else {
        (n / d).to_single()
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
