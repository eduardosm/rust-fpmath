use crate::traits::Float;

mod asin_acos;
mod atan;
mod cbrt;
mod div_pi;
mod exp;
mod exp10;
mod exp2;
mod gamma;
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

// GENERATE: consts f32
const PI: u32 = 0x40490FDB; // 3.1415927e0
const FRAC_PI_2: u32 = 0x3FC90FDB; // 1.5707964e0
const FRAC_PI_4: u32 = 0x3F490FDB; // 7.853982e-1
const FRAC_2_PI: u32 = 0x3F22F983; // 6.3661975e-1

impl<F: F32Like> crate::traits::FloatConsts<LikeF32> for F {
    #[inline]
    fn pi() -> Self {
        Self::from_raw(PI)
    }

    #[inline]
    fn frac_pi_2() -> Self {
        Self::from_raw(FRAC_PI_2)
    }

    #[inline]
    fn frac_pi_4() -> Self {
        Self::from_raw(FRAC_PI_4)
    }

    #[inline]
    fn frac_2_pi() -> Self {
        Self::from_raw(FRAC_2_PI)
    }
}
