use crate::traits::{Float, Int as _};

mod cbrt;
mod exp;
mod frexp;
mod gamma;
mod hyperbolic;
mod hypot;
mod inv_hyperbolic;
mod inv_trigonometric;
mod log;
mod pow;
mod round;
mod scalbn;
mod sqrt;
mod trigonometric;

pub(crate) use cbrt::{Cbrt, cbrt};
pub(crate) use exp::{Exp, exp, exp_m1, exp2, exp10};
pub(crate) use frexp::frexp;
pub(crate) use gamma::{Gamma, gamma, ln_gamma};
pub(crate) use hyperbolic::{Hyperbolic, cosh, sinh, sinh_cosh, tanh};
pub(crate) use hypot::{Hypot, hypot};
pub(crate) use inv_hyperbolic::{InvHyperbolic, acosh, asinh, atanh};
pub(crate) use inv_trigonometric::{
    InvTrigonometric, acos, acosd, acospi, asin, asind, asinpi, atan, atan2, atan2d, atan2pi,
    atand, atanpi,
};
pub(crate) use log::{Log, ln, ln_1p, log2, log10};
pub(crate) use pow::{Pow, pow, powi};
pub(crate) use round::{ceil, floor, round, trunc};
pub(crate) use scalbn::{scalbn, scalbn_medium};
pub(crate) use sqrt::{Sqrt, rsqrt_sqrt_32, sqrt};
pub(crate) use trigonometric::{
    Trigonometric, cos, cosd, cospi, sin, sin_cos, sind, sind_cosd, sinpi, sinpi_cospi, tan, tand,
    tanpi,
};

/// Returns whether `x` is positive, normal and finite.
#[inline]
fn is_pos_normal<F: Float>(x: F) -> bool {
    // sign and exponent
    let se = x.to_raw() >> F::MANT_BITS;
    se >= F::Raw::ONE && se < F::Raw::from(F::MAX_RAW_EXP)
}

/// Returns whether `x` is a non-zero integer, or infinity or NaN.
pub(crate) fn is_int<F: Float>(x: F) -> bool {
    let e = x.raw_exp();
    if e > F::EXP_OFFSET + F::RawExp::from(F::MANT_BITS) {
        true
    } else if e < F::EXP_OFFSET {
        false
    } else {
        let frac_shift = (F::EXP_OFFSET + F::RawExp::from(F::MANT_BITS)) - e;
        (x.to_raw() & !(F::Raw::MAX << frac_shift)) == F::Raw::ZERO
    }
}

fn is_odd_int<F: Float>(x: F) -> bool {
    let e = x.raw_exp();
    if e > F::EXP_OFFSET + F::RawExp::from(F::MANT_BITS) || e < F::EXP_OFFSET {
        // infinity, an even integer or only fractional part (less than 1)
        false
    } else {
        let frac_shift = (F::EXP_OFFSET + F::RawExp::from(F::MANT_BITS)) - e;
        if (x.to_raw() & !(F::Raw::MAX << frac_shift)) != F::Raw::ZERO {
            // not an integer
            false
        } else {
            ((x.mant() >> frac_shift) & F::Raw::ONE) == F::Raw::ONE
        }
    }
}

// like `is_odd_int`, but assumes that `x` is an integer
fn int_is_odd<F: Float>(x: F) -> bool {
    let e = x.raw_exp();
    if e > F::EXP_OFFSET + F::RawExp::from(F::MANT_BITS) {
        false
    } else {
        let frac_shift = (F::EXP_OFFSET + F::RawExp::from(F::MANT_BITS)) - e;
        ((x.mant() >> frac_shift) & F::Raw::ONE) == F::Raw::ONE
    }
}
