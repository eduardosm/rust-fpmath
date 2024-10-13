use crate::traits::{Float, Int as _};

mod acosh;
mod asin_acos;
mod asind_acosd;
mod asinh;
mod asinpi_acospi;
mod atan;
mod atand;
mod atanh;
mod atanpi;
mod cbrt;
mod ceil;
mod div_pi;
mod exp;
mod exp10;
mod exp2;
mod floor;
mod frexp;
mod gamma;
mod hypot;
mod log;
mod log10;
mod log2;
mod pow;
mod powi;
mod rad_to_deg;
mod reduce_90_deg;
mod reduce_half_mul_pi;
mod reduce_pi_2;
mod reduce_pi_2_large;
mod round;
mod scalbn;
mod sin_cos;
mod sind_cosd;
mod sinh_cosh;
mod sinpi_cospi;
mod sqrt;
mod tan;
mod tand;
mod tanh;
mod tanpi;
mod trunc;

pub(crate) use acosh::acosh;
pub(crate) use asin_acos::{acos, asin, AsinAcos};
pub(crate) use asind_acosd::{acosd, asind};
pub(crate) use asinh::asinh;
pub(crate) use asinpi_acospi::{acospi, asinpi};
pub(crate) use atan::{atan, atan2, Atan};
pub(crate) use atand::{atan2d, atand};
pub(crate) use atanh::atanh;
pub(crate) use atanpi::{atan2pi, atanpi};
pub(crate) use cbrt::{cbrt, Cbrt};
pub(crate) use ceil::ceil;
pub(crate) use div_pi::DivPi;
pub(crate) use exp::{exp, exp_m1, Exp};
pub(crate) use exp10::{exp10, Exp10};
pub(crate) use exp2::{exp2, Exp2};
pub(crate) use floor::floor;
pub(crate) use frexp::frexp;
pub(crate) use gamma::{lgamma, tgamma, Gamma};
pub(crate) use hypot::hypot;
pub(crate) use log::{log, log_1p, Log};
pub(crate) use log10::{log10, Log10};
pub(crate) use log2::{log2, Log2};
pub(crate) use pow::pow;
pub(crate) use powi::powi;
pub(crate) use rad_to_deg::RadToDeg;
pub(crate) use reduce_90_deg::{reduce_90_deg, Reduce90Deg};
pub(crate) use reduce_half_mul_pi::{reduce_half_mul_pi, ReduceHalfMulPi};
pub(crate) use reduce_pi_2::{reduce_pi_2, ReducePi2};
pub(crate) use round::{round, round_as_i_f};
pub(crate) use scalbn::{scalbn, scalbn_medium};
pub(crate) use sin_cos::{cos, sin, sin_cos, SinCos};
pub(crate) use sind_cosd::{cosd, sind, sind_cosd};
pub(crate) use sinh_cosh::{cosh, sinh, sinh_cosh, SinhCosh};
pub(crate) use sinpi_cospi::{cospi, sinpi, sinpi_cospi};
pub(crate) use sqrt::sqrt;
pub(crate) use tan::{tan, Tan};
pub(crate) use tand::tand;
pub(crate) use tanh::tanh;
pub(crate) use tanpi::tanpi;
pub(crate) use trunc::trunc;

fn is_int<F: Float>(x: F) -> bool {
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
