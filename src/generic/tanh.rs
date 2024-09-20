use super::exp::exp_split;
use super::sinh_cosh::{sinh_cosh_inner_common_1, sinh_cosh_inner_common_2};
use super::{Exp, SinhCosh};
use crate::traits::Int as _;

pub(crate) fn tanh<F: SinhCosh>(x: F) -> F {
    let e = x.raw_exp();
    if x >= F::expo2_hi_th() {
        // also handles x = inf
        F::one()
    } else if x <= -F::expo2_hi_th() {
        // also handles x = -inf
        -F::one()
    } else if e == F::MAX_RAW_EXP || e <= F::RawExp::ONE {
        // propagate NaN
        // or
        // very small, includes subnormal and zero
        // tanh(x) ~= x
        // also handles tanh(-0) = -0
        x
    } else {
        tanh_inner(x)
    }
}

fn tanh_inner<F: Exp>(x: F) -> F {
    // Split |x| into k, r_hi, r_lo such as:
    //  - |x| = k*ln(2) + r_hi + r_lo
    //  - k is an integer
    //  - |r_hi| <= 0.5*ln(2)
    let absx = x.abs();
    let (k, r_hi, r_lo) = exp_split(absx);

    if k > F::MANT_BITS.into() {
        F::one().copysign(x)
    } else {
        let (r_hi, r_lo) = F::norm_hi_lo_full(r_hi, r_lo);

        // t1a = exp(r_hi + r_lo) - 1 - r_hi
        // t1b = exp(-r_hi - r_lo) - 1 + r_hi
        let (t1a, t1b) = sinh_cosh_inner_common_1(r_hi, r_lo);

        // abss = |sinh(x)| = (exp(|x|) - exp(-|x|)) / 2
        // c = cosh(x) = (exp(|x|) + exp(-|x|)) / 2
        let (abss_hi, abss_lo, c_hi, c_lo) = sinh_cosh_inner_common_2(k, r_hi, t1a, t1b);

        // abst = |tanh(x)| = |sinh(x)| / cosh(x)

        let (n_hi, n_lo) = F::norm_hi_lo_splitted(abss_hi, abss_lo);
        let (d_hi, d_lo) = F::norm_hi_lo_splitted(c_hi, c_lo);

        let (q_hi, q_lo) = F::div_hi_lo(n_hi, n_lo, d_hi, d_lo);
        (q_hi + q_lo).copysign(x)
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::Float;
    use crate::FloatMath;

    fn test<F: Float + FloatMath>(hi_th: &str) {
        use crate::tanh;

        let test_nan = |arg: F| {
            let t = tanh(arg);
            assert_is_nan!(t);
        };

        let test_value = |arg: F, expected: F| {
            let t = tanh(arg);
            assert_total_eq!(t, expected);
        };

        let hi_th = F::parse(hi_th);

        test_nan(F::NAN);
        test_value(F::INFINITY, F::one());
        test_value(F::neg_infinity(), -F::one());
        test_value(hi_th, F::one());
        test_value(-hi_th, -F::one());
        test_value(hi_th + F::half(), F::one());
        test_value(-(hi_th + F::half()), -F::one());
        test_value(hi_th + F::one(), F::one());
        test_value(-(hi_th + F::one()), -F::one());
        test_value(F::ZERO, F::ZERO);
        test_value(-F::ZERO, -F::ZERO);
    }

    #[test]
    fn test_f32() {
        test::<f32>("89.5");
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f32() {
        test::<crate::SoftF32>("89.5");
    }

    #[test]
    fn test_f64() {
        test::<f64>("710.5");
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f64() {
        test::<crate::SoftF64>("710.5");
    }
}
