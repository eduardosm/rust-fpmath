use super::{F32Like, LikeF32};
use crate::double::SemiDouble;

// GENERATE: f32::log2::consts
const LOG2_E_HI: u32 = 0x3FB8A000; // 1.4423828e0
const LOG2_E_LO: u32 = 0x39A3B296; // 3.122284e-4

impl<F: F32Like> crate::generic::Log2<LikeF32> for F {
    #[inline]
    fn log2_e_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(Self::from_raw(LOG2_E_HI), Self::from_raw(LOG2_E_LO))
    }
}
