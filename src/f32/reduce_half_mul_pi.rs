use super::{F32Like, LikeF32};
use crate::double::SemiDouble;

// GENERATE: reduce_half_mul_pi::consts f32
const PI_HI: u32 = 0x40490000; // 3.140625e0
const PI_LO: u32 = 0x3A7DAA22; // 9.676536e-4

impl<F: F32Like> crate::generic::ReduceHalfMulPi<LikeF32> for F {
    #[inline]
    fn pi_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(Self::from_raw(PI_HI), Self::from_raw(PI_LO))
    }
}
