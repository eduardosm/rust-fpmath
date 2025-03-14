use super::{F32Like, LikeF32};

// GENERATE: exp2::consts f32
const LN_2: u32 = 0x3F317218; // 6.931472e-1

impl<F: F32Like> crate::generic::Exp2<LikeF32> for F {
    #[inline]
    fn ln_2() -> Self {
        Self::from_raw(LN_2)
    }

    #[inline]
    fn exp2_lo_th() -> Self {
        Self::cast_from(-151i16)
    }

    #[inline]
    fn exp2_hi_th() -> Self {
        Self::cast_from(129i16)
    }
}
