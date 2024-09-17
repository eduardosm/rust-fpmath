use super::{F32Like, LikeF32};

// Generated with `./run-generator.sh f32::exp2::consts`
const LN_2: u32 = 0x3F317218; // 6.931472e-1
const LN_2_HI: u32 = 0x3F317000; // 6.9311523e-1
const LN_2_LO: u32 = 0x3805FDF4; // 3.1946183e-5

impl<F: F32Like> crate::generic::Exp2<LikeF32> for F {
    #[inline]
    fn ln_2() -> Self {
        Self::from_raw(LN_2)
    }

    #[inline]
    fn ln_2_hi() -> Self {
        Self::from_raw(LN_2_HI)
    }

    #[inline]
    fn ln_2_lo() -> Self {
        Self::from_raw(LN_2_LO)
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
