use crate::double::{DenormDouble, SemiDouble};
use crate::traits::{Float, Int as _};

pub(crate) fn sqrt<F: Float>(x: F) -> F {
    let (y, edelta) = x.normalize_arg();
    let yexp = y.raw_exp();
    if yexp == F::RawExp::ZERO {
        // sqrt(±0) = ±0
        y
    } else if y.sign() {
        // x < 0, sqrt(x) = NaN
        F::NAN
    } else if yexp == F::MAX_RAW_EXP {
        // propagate infinity or NaN
        y
    } else {
        sqrt_inner(y, edelta)
    }
}

fn sqrt_inner<F: Float>(x: F, edelta: F::Exp) -> F {
    // Split x * 2^edelta = 2^(k - MANT_BITS) * m such as
    // * k is an even integer
    // * m is a positive integer
    // * 1 <= m * 2^(-MANT_BITS) < 4
    let (k, m) = sqrt_split(x, edelta);

    // Calculate sqrt(m) bit-by-bit with one extra bit
    let mut rem = m << 1; // extra bit
    let mut res = F::Raw::ZERO;
    let mut s = F::Raw::ZERO;
    let mut d = F::Raw::ONE << (F::MANT_BITS + 1);

    while d != F::Raw::ZERO {
        let t = s + d;
        if rem >= t {
            s += d + d;
            rem -= t;
            res += d;
        }
        rem <<= 1;
        d >>= 1;
    }

    // Round to nearest and remove extra bit
    res = (res + F::Raw::ONE) >> 1;

    // Build result with k/2 as exponent and res as mantissa
    let raw_exp = F::exp_to_raw_exp(k >> 1) - F::RawExp::ONE;
    let yraw = (F::Raw::from(raw_exp) << F::MANT_BITS) + res;
    F::from_raw(yraw)
}

fn sqrt_split<F: Float>(x: F, edelta: F::Exp) -> (F::Exp, F::Raw) {
    // Split x * 2^edelta = 2^(k - MANT_BITS) * m

    let k = x.exponent() + edelta;
    let m = x.mant();

    if (k & F::Exp::ONE) == F::Exp::ZERO {
        // exponent is even, 1 <= m * 2^(-MANT_BITS) < 2
        (k, m)
    } else {
        // exponent is odd, 2 <= m * 2^(-MANT_BITS) < 4
        (k - F::Exp::ONE, m << 1)
    }
}

/// Calculates `2 * sqrt(x)` with extended precision
pub(super) fn two_hi_lo_sqrt_inner<F: Float>(x: F) -> DenormDouble<F> {
    // Improve accuracy with a single Newton iteration
    // y = sqrt(x)
    // sqrt(x)_hi + sqrt(x)_lo = (y * y + x) / (2 * y)

    let y = sqrt_inner(x, F::Exp::ZERO);
    let y = SemiDouble::new(y);

    let y2 = y.square();

    SemiDouble::new_qadd12(x, y2) / y
}

/// Calculates `sqrt(x_hi + x_lo)` with extended precision
pub(super) fn hi_lo_sqrt_hi_lo_inner<F: Float>(x: DenormDouble<F>) -> DenormDouble<F> {
    // Improve accuracy with a single Newton iteration
    // y = sqrt(x)
    // sqrt(x)_hi + sqrt(x)_lo = (y * y + x_hi + x_lo) / (2 * y)

    let y = sqrt_inner(x.hi(), F::Exp::ZERO);
    let y = SemiDouble::new(y);

    let y2 = y.square();

    SemiDouble::new_qadd22(x, y2) / y.pmul1(F::two())
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
        use crate::{scalbn, sqrt};

        assert_is_nan!(sqrt(F::NAN));
        assert_is_nan!(sqrt(F::neg_infinity()));
        assert_is_nan!(sqrt(-F::one()));
        assert_total_eq!(sqrt(F::INFINITY), F::INFINITY);
        assert_total_eq!(sqrt(F::ZERO), F::ZERO);
        assert_total_eq!(sqrt(-F::ZERO), -F::ZERO);

        let min_normal_exp: i32 = F::MIN_NORMAL_EXP.into();
        let e_limit = (-min_normal_exp) / 2;
        for e in (-e_limit)..e_limit {
            for m in full_e_mants.clone() {
                let x = scalbn(F::cast_from(m), e - m.ilog2() as i32);
                let x2 = x * x;
                assert_total_eq!(sqrt(x2), x);
            }
        }
        for e in extra_e {
            for m in extra_e_mants.clone() {
                let x = scalbn(F::cast_from(m), e - m.ilog2() as i32);
                let x2 = x * x;
                assert_total_eq!(sqrt(x2), x);
            }
        }
    }

    #[test]
    fn test_f32() {
        test::<f32>(0x800..=0xFFF, core::iter::empty(), core::iter::empty());
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f32() {
        test::<crate::SoftF32>(0x800..=0xFFF, core::iter::empty(), core::iter::empty());
    }

    #[test]
    fn test_f64() {
        test::<f64>(0x800..=0xFFF, [-511, 0, 511].into_iter(), 0x8000..=0xFFFF);
    }

    #[cfg(feature = "soft-float")]
    #[test]
    fn test_soft_f64() {
        test::<crate::SoftF64>(0x800..=0xFFF, [-511, 0, 511].into_iter(), 0x8000..=0xFFFF);
    }
}
