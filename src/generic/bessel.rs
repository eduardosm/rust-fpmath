use crate::traits::{Float, Int as _};

pub(crate) trait Bessel: Float {
    fn j0_finite(x: Self) -> Self;

    fn j1_finite(x: Self) -> Self;

    fn y0_finite(x: Self) -> Self;

    fn y1_finite(x: Self) -> Self;

    /// `jn(n, x)` for `n >= 2` and finite `x > 0`.
    fn jn_finite(n: u32, x: Self) -> Self;

    /// `yn(n, x)` for `n >= 2` and finite `x > 0`.
    fn yn_finite(n: u32, x: Self) -> Self;
}

pub(crate) fn j0<F: Bessel>(x: F) -> F {
    if x.raw_exp() == F::MAX_RAW_EXP {
        if x.raw_mant() == F::Raw::ZERO {
            // j0(±inf) = 0
            F::ZERO
        } else {
            // propagate NaN
            x
        }
    } else {
        // also handles zero and subnormals, j0(x) ~= 1
        F::j0_finite(x)
    }
}

pub(crate) fn j1<F: Bessel>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // j1(±0) = ±0
        x
    } else if e == F::MAX_RAW_EXP {
        if x.raw_mant() == F::Raw::ZERO {
            // j1(±inf) = ±0
            F::ZERO.copysign(x)
        } else {
            // propagate NaN
            x
        }
    } else {
        // also handles subnormals, j1(x) ~= x / 2
        F::j1_finite(x)
    }
}

pub(crate) fn y0<F: Bessel>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // y0(±0) = -inf
        F::NEG_INFINITY
    } else if e == F::MAX_RAW_EXP && x.raw_mant() != F::Raw::ZERO {
        // propagate NaN
        x
    } else if x.is_sign_negative() {
        // x < 0
        F::NAN
    } else if e == F::MAX_RAW_EXP {
        // y0(inf) = 0
        F::ZERO
    } else {
        // also handles subnormals
        F::y0_finite(x)
    }
}

pub(crate) fn y1<F: Bessel>(x: F) -> F {
    let e = x.raw_exp();
    if e == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO {
        // y1(±0) = -inf
        F::NEG_INFINITY
    } else if e == F::MAX_RAW_EXP && x.raw_mant() != F::Raw::ZERO {
        // propagate NaN
        x
    } else if x.is_sign_negative() {
        // x < 0
        F::NAN
    } else if e == F::MAX_RAW_EXP {
        // y1(inf) = 0
        F::ZERO
    } else {
        // also handles subnormals
        F::y1_finite(x)
    }
}

pub(crate) fn jn<F: Bessel>(n: i32, x: F) -> F {
    let e = x.raw_exp();
    if e == F::MAX_RAW_EXP && x.raw_mant() != F::Raw::ZERO {
        // propagate NaN
        return x;
    }
    if n == 0 {
        return j0(x);
    }

    // jn(-n, x) = (-1)^n * jn(n, x) and jn(n, -x) = (-1)^n * jn(n, x), so it
    // is calculated with the absolute values, and then the sign is applied
    let order = n.unsigned_abs();
    let r = if order == 1 {
        j1(x.abs())
    } else if e == F::MAX_RAW_EXP || (e == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO) {
        // jn(n, 0) = 0 and jn(n, inf) = 0
        F::ZERO
    } else {
        // also handles subnormals
        F::jn_finite(order, x.abs())
    };
    let negate = (order & 1) != 0 && ((n < 0) != x.is_sign_negative());
    if negate { -r } else { r }
}

pub(crate) fn yn<F: Bessel>(n: i32, x: F) -> F {
    let e = x.raw_exp();
    let is_zero = e == F::RawExp::ZERO && x.raw_mant() == F::Raw::ZERO;
    if e == F::MAX_RAW_EXP && x.raw_mant() != F::Raw::ZERO {
        // propagate NaN
        return x;
    }
    if n == 0 {
        return y0(x);
    }
    if x.is_sign_negative() && !is_zero {
        // x < 0
        return F::NAN;
    }

    // yn(-n, x) = (-1)^n * yn(n, x), so it is calculated with the absolute
    // value of `n`, and then the sign is applied
    let order = n.unsigned_abs();
    let r = if order == 1 {
        y1(x)
    } else if is_zero {
        // yn(n, ±0) = -inf
        F::NEG_INFINITY
    } else if e == F::MAX_RAW_EXP {
        // yn(n, inf) = 0
        F::ZERO
    } else {
        // also handles subnormals
        F::yn_finite(order, x)
    };
    let negate = (order & 1) != 0 && n < 0;
    if negate { -r } else { r }
}

#[cfg(test)]
mod tests {
    use crate::FloatMath;
    use crate::traits::{Float, Int as _};

    fn test_j0<F: Float + FloatMath>() {
        use crate::j0;

        assert_is_nan!(j0(F::NAN));
        assert_total_eq!(j0(F::INFINITY), F::ZERO);
        assert_total_eq!(j0(F::NEG_INFINITY), F::ZERO);
        assert_total_eq!(j0(F::ZERO), F::ONE);
        assert_total_eq!(j0(-F::ZERO), F::ONE);
    }

    fn test_j1<F: Float + FloatMath>() {
        use crate::j1;

        assert_is_nan!(j1(F::NAN));
        assert_total_eq!(j1(F::INFINITY), F::ZERO);
        assert_total_eq!(j1(F::NEG_INFINITY), -F::ZERO);
        assert_total_eq!(j1(F::ZERO), F::ZERO);
        assert_total_eq!(j1(-F::ZERO), -F::ZERO);
    }

    fn test_y0_y1<F: Float + FloatMath>(y: fn(F) -> F) {
        assert_is_nan!(y(F::NAN));
        assert_total_eq!(y(F::INFINITY), F::ZERO);
        assert_is_nan!(y(F::NEG_INFINITY));
        assert_total_eq!(y(F::ZERO), F::NEG_INFINITY);
        assert_total_eq!(y(-F::ZERO), F::NEG_INFINITY);
        assert_is_nan!(y(-F::ONE));
        assert_is_nan!(y(-F::LARGEST));
        assert_is_nan!(y(-F::from_raw(F::Raw::ONE)));
    }

    fn test_jn<F: Float + FloatMath>() {
        use crate::{j0, j1, jn};

        let f = F::parse;

        for n in [i32::MIN, -1000, -3, -2, -1, 0, 1, 2, 3, 1000, i32::MAX] {
            assert_is_nan!(jn(n, F::NAN));

            // The results are zero, with the sign of `x` when `n` is odd and
            // positive
            let sign = if n % 2 == 0 {
                [F::ZERO, F::ZERO]
            } else if n > 0 {
                [F::ZERO, -F::ZERO]
            } else {
                [-F::ZERO, F::ZERO]
            };
            assert_total_eq!(jn(n, F::INFINITY), sign[0]);
            assert_total_eq!(jn(n, F::NEG_INFINITY), sign[1]);
            if n != 0 {
                assert_total_eq!(jn(n, F::ZERO), sign[0]);
                assert_total_eq!(jn(n, -F::ZERO), sign[1]);
            }
        }
        assert_total_eq!(jn(0, F::ZERO), F::ONE);
        assert_total_eq!(jn(0, -F::ZERO), F::ONE);

        for x in ["0.25", "1", "2.5", "7", "30", "1000", "1e30"].map(f) {
            // Orders 0 and 1
            assert_total_eq!(jn(0, x), j0(x));
            assert_total_eq!(jn(0, -x), j0(x));
            assert_total_eq!(jn(1, x), j1(x));
            assert_total_eq!(jn(1, -x), -j1(x));
            assert_total_eq!(jn(-1, x), -j1(x));
            assert_total_eq!(jn(-1, -x), j1(x));

            // Symmetry
            for n in [2, 3, 10, 21] {
                let r = jn(n, x);
                let odd = if n % 2 == 0 { r } else { -r };
                assert_total_eq!(jn(n, -x), odd);
                assert_total_eq!(jn(-n, x), odd);
                assert_total_eq!(jn(-n, -x), r);
            }
        }

        // The results underflow when the order is much larger than `x`
        for n in [1000, 1001, i32::MAX] {
            assert_total_eq!(jn(n, F::ONE), F::ZERO);
            assert_total_eq!(jn(n, -F::ONE), if n % 2 == 0 { F::ZERO } else { -F::ZERO });
            assert_total_eq!(jn(-n, F::ONE), if n % 2 == 0 { F::ZERO } else { -F::ZERO });
        }
        assert_total_eq!(jn(i32::MIN, F::ONE), F::ZERO);
        assert_total_eq!(jn(2, F::from_raw(F::Raw::ONE)), F::ZERO);
    }

    fn test_yn<F: Float + FloatMath>() {
        use crate::{y0, y1, yn};

        let f = F::parse;

        for n in [i32::MIN, -1000, -3, -2, -1, 0, 1, 2, 3, 1000, i32::MAX] {
            // The results have the opposite sign when `n` is odd and negative
            let negate = n < 0 && n % 2 != 0;
            let sign = |v: F| if negate { -v } else { v };

            assert_is_nan!(yn(n, F::NAN));
            assert_total_eq!(yn(n, F::INFINITY), sign(F::ZERO));
            assert_is_nan!(yn(n, F::NEG_INFINITY));
            assert_total_eq!(yn(n, F::ZERO), sign(F::NEG_INFINITY));
            assert_total_eq!(yn(n, -F::ZERO), sign(F::NEG_INFINITY));
            assert_is_nan!(yn(n, -F::ONE));
            assert_is_nan!(yn(n, -F::LARGEST));
            assert_is_nan!(yn(n, -F::from_raw(F::Raw::ONE)));

            // The results overflow when the order is much larger than `x`
            // (and with tiny `x`)
            if n.unsigned_abs() >= 1000 {
                assert_total_eq!(yn(n, F::ONE), sign(F::NEG_INFINITY));
            }
            if n.unsigned_abs() >= 2 {
                assert_total_eq!(yn(n, F::from_raw(F::Raw::ONE)), sign(F::NEG_INFINITY));
            }
        }

        for x in ["0.25", "1", "2.5", "7", "30", "1000", "1e30"].map(f) {
            // Orders 0 and 1
            assert_total_eq!(yn(0, x), y0(x));
            assert_total_eq!(yn(1, x), y1(x));
            assert_total_eq!(yn(-1, x), -y1(x));

            // Symmetry
            for n in [2, 3, 10, 21] {
                let r = yn(n, x);
                assert_total_eq!(yn(-n, x), if n % 2 == 0 { r } else { -r });
            }
        }
    }

    #[test]
    fn test_f32() {
        test_j0::<f32>();
        test_j1::<f32>();
        test_y0_y1::<f32>(crate::y0);
        test_y0_y1::<f32>(crate::y1);
        test_jn::<f32>();
        test_yn::<f32>();
    }

    #[test]
    fn test_f64() {
        test_j0::<f64>();
        test_j1::<f64>();
        test_y0_y1::<f64>(crate::y0);
        test_y0_y1::<f64>(crate::y1);
        test_jn::<f64>();
        test_yn::<f64>();
    }
}
