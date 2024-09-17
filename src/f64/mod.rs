use crate::traits::Float;

mod asin_acos;
mod atan;
mod cbrt;
mod div_pi;
mod exp;
mod exp10;
mod exp2;
mod log;
mod log10;
mod log2;
mod rad_to_deg;
mod reduce_90_deg;
mod reduce_half_mul_pi;
mod reduce_pi_2;
mod sin_cos;
mod sinh_cosh;
mod tan;

pub(crate) struct LikeF64;

pub(crate) trait F64Like: Float<Like = LikeF64, Raw = u64, RawExp = u16, Exp = i16> {}
