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

// Generated with `./run-generator.sh f64::consts`
const PI: u64 = 0x400921FB54442D18; // 3.141592653589793e0
const FRAC_PI_2: u64 = 0x3FF921FB54442D18; // 1.5707963267948966e0
const FRAC_PI_4: u64 = 0x3FE921FB54442D18; // 7.853981633974483e-1
const FRAC_1_PI: u64 = 0x3FD45F306DC9C883; // 3.183098861837907e-1
const FRAC_2_PI: u64 = 0x3FE45F306DC9C883; // 6.366197723675814e-1

impl<F: F64Like> crate::traits::FloatConsts<LikeF64> for F {
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
    fn frac_1_pi() -> Self {
        Self::from_raw(FRAC_1_PI)
    }

    #[inline]
    fn frac_2_pi() -> Self {
        Self::from_raw(FRAC_2_PI)
    }
}
