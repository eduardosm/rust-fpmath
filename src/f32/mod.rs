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

pub(crate) struct LikeF32;

pub(crate) trait F32Like: Float<Like = LikeF32, Raw = u32, RawExp = u8, Exp = i16> {}
