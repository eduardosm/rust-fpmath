use crate::double::{DenormDouble, NormDouble, SemiDouble};
use crate::generic::{reduce_pi_2, ReducePi2};
use crate::traits::{Float, Int as _, Like};

pub(crate) trait SinCos<L = Like<Self>>: Float {
    fn frac_1_6_ex() -> SemiDouble<Self>;

    /// Calculates `(sin(x) - x - x^3 * K3, K3)`
    ///
    /// Where:
    /// * `x2 = x^2`
    /// * `x5 = x^5`
    /// * `K3 ~= -1/6`
    fn sin_poly(x2: Self, x5: Self) -> (Self, Self);

    /// Calculates `sin(x) - x + x^3 * 1/6`
    ///
    /// Where:
    /// * `x2 = x^2`
    /// * `x5 = x^5`
    fn sin_poly_ex(x2: Self, x5: Self) -> Self;

    /// Calculates `cos(x) + 0.5 * x^2 - 1`
    ///
    /// Where:
    /// * `x2 = x^2`
    /// * `x4 = x^4`
    fn cos_poly(x2: Self, x4: Self) -> Self;
}

pub(crate) fn sin<F: SinCos + ReducePi2>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        // sin(inf or nan) = nan
        F::NAN
    } else if e == F::RawExp::ZERO {
        // subnormal or zero, sin(x) ~= x
        // also handles sin(-0) = -0
        x
    } else {
        let (n, y_hi, y_lo) = reduce_pi_2(x);

        match n {
            0 => sin_inner(y_hi, y_lo),
            1 => cos_inner(y_hi, y_lo),
            2 => -sin_inner(y_hi, y_lo),
            3 => -cos_inner(y_hi, y_lo),
            _ => unreachable!(),
        }
    }
}

pub(crate) fn cos<F: SinCos + ReducePi2>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        // cos(inf or nan) = nan
        F::NAN
    } else if e == F::RawExp::ZERO {
        // subnormal or zero, cos(x) ~= 1
        F::one()
    } else {
        let (n, y_hi, y_lo) = reduce_pi_2(x);

        match n {
            0 => cos_inner(y_hi, y_lo),
            1 => -sin_inner(y_hi, y_lo),
            2 => -cos_inner(y_hi, y_lo),
            3 => sin_inner(y_hi, y_lo),
            _ => unreachable!(),
        }
    }
}

pub(crate) fn sin_cos<F: SinCos + ReducePi2>(x: F) -> (F, F) {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        // sin(inf or nan) = nan
        // cos(inf or nan) = nan
        (F::NAN, F::NAN)
    } else if e == F::RawExp::ZERO {
        // subnormal or zero
        // sin(x) ~= x
        // cos(x) ~= 1
        // also handles sin(-0) = -0
        (x, F::one())
    } else {
        let (n, y_hi, y_lo) = reduce_pi_2(x);

        let sin = sin_inner(y_hi, y_lo);
        let cos = cos_inner(y_hi, y_lo);
        match n {
            0 => (sin, cos),
            1 => (cos, -sin),
            2 => (-sin, -cos),
            3 => (-cos, sin),
            _ => unreachable!(),
        }
    }
}

/// Calculates `sin(x_hi + x_lo)`, where
/// `x_lo` is very small and `|x_hi| <= π/4`
pub(super) fn sin_inner<F: SinCos>(x_hi: F, x_lo: F) -> F {
    // sin(x_hi + x_lo) = sin(x_hi) * cos(x_lo) + cos(x_hi) * sin(x_lo)
    // x_lo is small, so sin(x_lo) ~= x_lo and cos(x_lo) ~= 1,
    // then sin(x_hi + x_lo) ~= sin(x_hi) + cos(x_hi) * x_lo
    //
    // sin(x) is calculated with a polynomial.
    // cos(x) * x_lo ~= (1 - 0.5*x^2) * x_lo

    let x2 = x_hi * x_hi;
    let x3 = x2 * x_hi;
    let x5 = x3 * x2;

    // t1 = sin(x_hi) - x_hi - x_hi^3 * k3
    let (t1, k3) = F::sin_poly(x2, x5);

    // sin(x_hi + x_lo) ~= sin(x) + (1 - 0.5 * x^2) * x_lo
    //                  = x + t1 + x^3 * k3 + (1 - 0.5 * x^2) * x_lo
    x_hi + (x3 * k3 + (t1 + (x_lo - F::half() * x2 * x_lo)))
}

pub(super) fn hi_lo_sin_inner<F: SinCos>(x: NormDouble<F>) -> DenormDouble<F> {
    // sin(x) is calculated with a polynomial.

    let x_semi = x.to_semi();
    let x2 = x_semi.square();
    let x2_single = x2.to_single();
    let x3 = x2.to_semi() * x_semi;
    let x5 = x3.to_single() * x2_single;

    // t1 = sin(x) - x + x^3 / 6
    let t1 = F::sin_poly_ex(x2_single, x5);

    // sin(x) = t1 + x - x^3 / 6
    let x3k3 = x3.to_semi() * (-F::frac_1_6_ex());
    x.to_denorm().qadd2(x3k3 + t1)
}

/// Calculates `cos(x_hi + x_lo)`, where
/// `x_lo` is very small and `|x_hi| <= π/4`
pub(super) fn cos_inner<F: SinCos>(x_hi: F, x_lo: F) -> F {
    // cos(x_hi + x_lo) = cos(x_hi) * cos(x_lo) - sin(x_hi) * sin(x_lo)
    // x_lo is small, so sin(x_lo) ~= x_lo and cos(x_lo) ~= 1
    // then
    // cos(x_hi + x_lo) ~= cos(x_hi) - sin(x_hi) * x_lo
    //
    // cos(x_hi) is calculated with a polynomial.
    // sin(x_hi) * x_lo ~= x_hi * x_lo

    // t1 = cos(x_hi) + 0.5 * x_hi^2 - 1
    let x2 = x_hi * x_hi;
    let x4 = x2 * x2;
    let t1 = F::cos_poly(x2, x4);

    // cos(x_hi + x_lo) = t1 + 1 - 0.5 * x^2 - x_hi * x_lo
    let t2 = DenormDouble::new_qsub11(F::one(), F::half() * x2);
    t2.lsub(x_hi * x_lo).ladd(t1).to_single()
}

pub(super) fn hi_lo_cos_inner<F: SinCos>(x: NormDouble<F>) -> DenormDouble<F> {
    // cos(x) is calculated with a polynomial.

    // t1 = cos(x) + 0.5 * x^2 - 1
    let x2 = x.to_semi().square();
    let x2_single = x2.to_single();
    let x4 = x2_single * x2_single;
    let t1 = F::cos_poly(x2_single, x4);

    // t2 = 1 - 0.5 * x^2
    let t2 = DenormDouble::new_qsub12(F::one(), x2.pmul1(F::half()));

    // cos(x) = t1 + t2
    t2.qadd1(t1)
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test<F: Float + FloatMath>() {
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
        test_nan(F::neg_infinity());
        test_value(F::ZERO, F::ZERO, F::one());
        test_value(-F::ZERO, -F::ZERO, F::one());
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
