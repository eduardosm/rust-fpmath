use super::{F32Like, LikeF32};

// Generated with `./run-generator.sh f32::exp10::consts`
const LOG2_10: u32 = 0x40549A78; // 3.321928e0
const LOG10_2_HI: u32 = 0x3E9A2000; // 3.010254e-1
const LOG10_2_LO: u32 = 0x369A84FC; // 4.605039e-6
const LN_10: u32 = 0x40135D8E; // 2.3025851e0
const LN_10_HI: u32 = 0x40135000; // 2.3017578e0
const LN_10_LO: u32 = 0x3A58DDDB; // 8.272805e-4

impl<F: F32Like> crate::generic::Exp10<LikeF32> for F {
    #[inline]
    fn log2_10() -> Self {
        Self::from_raw(LOG2_10)
    }

    #[inline]
    fn log10_2_hi() -> Self {
        Self::from_raw(LOG10_2_HI)
    }

    #[inline]
    fn log10_2_lo() -> Self {
        Self::from_raw(LOG10_2_LO)
    }

    #[inline]
    fn ln_10() -> Self {
        Self::from_raw(LN_10)
    }

    #[inline]
    fn ln_10_hi() -> Self {
        Self::from_raw(LN_10_HI)
    }

    #[inline]
    fn ln_10_lo() -> Self {
        Self::from_raw(LN_10_LO)
    }

    #[inline]
    fn exp10_lo_th() -> Self {
        Self::cast_from(-46i16)
    }

    #[inline]
    fn exp10_hi_th() -> Self {
        Self::cast_from(39i16)
    }
}
