use super::{F64Like, LikeF64};
use crate::double::SemiDouble;

// GENERATE: reduce_half_mul_pi::consts f64
const PI_HI: u64 = 0x400921FB50000000; // 3.1415926218032837e0
const PI_LO: u64 = 0x3E6110B4611A6263; // 3.178650954705639e-8

impl<F: F64Like> crate::generic::ReduceHalfMulPi<LikeF64> for F {
    #[inline]
    fn pi_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(Self::from_raw(PI_HI), Self::from_raw(PI_LO))
    }
}
