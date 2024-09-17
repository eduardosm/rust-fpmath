use super::{F32Like, LikeF32};

// Generated with `./run-generator.sh f32::log2::consts`
const LOG2_E_HI: u32 = 0x3FB8A000; // 1.4423828e0
const LOG2_E_LO: u32 = 0x39A3B296; // 3.122284e-4

impl<F: F32Like> crate::generic::Log2<LikeF32> for F {
    #[inline]
    fn log2_e_hi() -> Self {
        Self::from_raw(LOG2_E_HI)
    }

    #[inline]
    fn log2_e_lo() -> Self {
        Self::from_raw(LOG2_E_LO)
    }
}
