use super::{F32Like, LikeF32};
use crate::double::SemiDouble;

// GENERATE: log10::consts f32
const LOG10_E_HI: u32 = 0x3EDE5000; // 4.342041e-1
const LOG10_E_LO: u32 = 0x38BD8A93; // 9.038034e-5
const LOG10_2_HI: u32 = 0x3E9A2000; // 3.010254e-1
const LOG10_2_LO: u32 = 0x369A84FC; // 4.605039e-6

impl<F: F32Like> crate::generic::Log10<LikeF32> for F {
    #[inline]
    fn log10_e_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(Self::from_raw(LOG10_E_HI), Self::from_raw(LOG10_E_LO))
    }

    #[inline]
    fn log10_2_hi() -> Self {
        Self::from_raw(LOG10_2_HI)
    }

    #[inline]
    fn log10_2_lo() -> Self {
        Self::from_raw(LOG10_2_LO)
    }
}
