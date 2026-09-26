use crate::traits::{Float, Int as _};

pub(crate) trait Hyperbolic: Float {
    fn sinh_finite(x: Self) -> Self;

    fn cosh_finite(x: Self) -> Self;

    fn sinh_cosh_finite(x: Self) -> (Self, Self);

    fn tanh_finite(x: Self) -> Self;
}

pub(crate) fn sinh<F: Hyperbolic>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP || e <= F::RawExp::ONE {
        // propagate NaN
        // or
        // sinh(±inf) = ±inf
        // or
        // very small, includes subnormal and zero
        // sinh(x) ~= x
        // also handles sinh(-0) = -0
        x
    } else {
        F::sinh_finite(x)
    }
}

pub(crate) fn cosh<F: Hyperbolic>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        if x.raw_mant() == F::Raw::ZERO {
            // cosh(±inf) = +inf
            F::INFINITY
        } else {
            // Propagate NaN
            x
        }
    } else if e <= F::RawExp::ONE {
        // very small, includes subnormal and zero
        // cosh(x) ~= 1
        F::ONE
    } else {
        F::cosh_finite(x)
    }
}

pub(crate) fn sinh_cosh<F: Hyperbolic>(x: F) -> (F, F) {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        if x.raw_mant() == F::Raw::ZERO {
            // sinh(±inf) = ±inf
            // cosh(±inf) = +inf
            (x, F::INFINITY)
        } else {
            // Propagate NaN
            (x, x)
        }
    } else if e <= F::RawExp::ONE {
        // very small, includes subnormal and zero
        // sinh(x) ~= x
        // cosh(x) ~= 1
        // also handles sinh(-0) = -0
        (x, F::ONE)
    } else {
        F::sinh_cosh_finite(x)
    }
}

pub(crate) fn tanh<F: Hyperbolic>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP {
        if x.raw_mant() == F::Raw::ZERO {
            // tanh(±inf) = ±1
            F::ONE.copysign(x)
        } else {
            // propagate NaN
            x
        }
    } else if e <= F::RawExp::ONE {
        // very small, includes subnormal and zero
        // tanh(x) ~= x
        // also handles tanh(-0) = -0
        x
    } else {
        F::tanh_finite(x)
    }
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::Float;

    fn test_sinh_cosh<F: Float + FloatMath>(hi_th: F) {
        use crate::{cosh, sinh, sinh_cosh};

        let test_nan = |arg: F| {
            let sinh1 = sinh(arg);
            let cosh1 = cosh(arg);
            let (sinh2, cosh2) = sinh_cosh(arg);
            assert_is_nan!(sinh1);
            assert_is_nan!(cosh1);
            assert_is_nan!(sinh2);
            assert_is_nan!(cosh2);
        };

        let test_value = |arg: F, expected_sinh: F, expected_cosh: F| {
            let sinh1 = sinh(arg);
            let cosh1 = cosh(arg);
            let (sinh2, cosh2) = sinh_cosh(arg);
            assert_total_eq!(sinh1, expected_sinh);
            assert_total_eq!(cosh1, expected_cosh);
            assert_total_eq!(sinh2, expected_sinh);
            assert_total_eq!(cosh2, expected_cosh);
        };

        test_nan(F::NAN);
        test_value(F::INFINITY, F::INFINITY, F::INFINITY);
        test_value(F::NEG_INFINITY, F::NEG_INFINITY, F::INFINITY);
        test_value(hi_th, F::INFINITY, F::INFINITY);
        test_value(-hi_th, F::NEG_INFINITY, F::INFINITY);
        test_value(hi_th + F::HALF, F::INFINITY, F::INFINITY);
        test_value(-(hi_th + F::HALF), F::NEG_INFINITY, F::INFINITY);
        test_value(hi_th + F::ONE, F::INFINITY, F::INFINITY);
        test_value(-(hi_th + F::ONE), F::NEG_INFINITY, F::INFINITY);
        test_value(F::ZERO, F::ZERO, F::ONE);
        test_value(-F::ZERO, -F::ZERO, F::ONE);
    }

    fn test_tanh<F: Float + FloatMath>(hi_th: F) {
        use crate::tanh;

        let test_nan = |arg: F| {
            let t = tanh(arg);
            assert_is_nan!(t);
        };

        let test_value = |arg: F, expected: F| {
            let t = tanh(arg);
            assert_total_eq!(t, expected);
        };

        test_nan(F::NAN);
        test_value(F::INFINITY, F::ONE);
        test_value(F::NEG_INFINITY, -F::ONE);
        test_value(hi_th, F::ONE);
        test_value(-hi_th, -F::ONE);
        test_value(hi_th + F::HALF, F::ONE);
        test_value(-(hi_th + F::HALF), -F::ONE);
        test_value(hi_th + F::ONE, F::ONE);
        test_value(-(hi_th + F::ONE), -F::ONE);
        test_value(F::ZERO, F::ZERO);
        test_value(-F::ZERO, -F::ZERO);
    }

    #[test]
    fn test_f32() {
        test_sinh_cosh::<f32>(89.5);
        test_tanh::<f32>(89.5);
    }

    #[test]
    fn test_f64() {
        test_sinh_cosh::<f64>(710.5);
        test_tanh::<f64>(710.5);
    }
}
