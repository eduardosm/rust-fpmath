// GENERATE: exp10::consts f32
const LOG2_10: u32 = 0x40549A78; // 3.321928e0
const LOG10_2_HI: u32 = 0x3E9A2000; // 3.010254e-1
const LOG10_2_LO: u32 = 0x369A84FC; // 4.605039e-6
const LN_10: u32 = 0x40135D8E; // 2.3025851e0
const LN_10_HI: u32 = 0x40135000; // 2.3017578e0
const LN_10_LO: u32 = 0x3A58DDDB; // 8.272805e-4

impl crate::generic::Exp10 for f32 {
    #[inline]
    fn log2_10() -> Self {
        f32::from_bits(LOG2_10)
    }

    #[inline]
    fn log10_2_hi() -> Self {
        f32::from_bits(LOG10_2_HI)
    }

    #[inline]
    fn log10_2_lo() -> Self {
        f32::from_bits(LOG10_2_LO)
    }

    #[inline]
    fn ln_10() -> Self {
        f32::from_bits(LN_10)
    }

    #[inline]
    fn ln_10_hi() -> Self {
        f32::from_bits(LN_10_HI)
    }

    #[inline]
    fn ln_10_lo() -> Self {
        f32::from_bits(LN_10_LO)
    }

    #[inline]
    fn exp10_lo_th() -> Self {
        -46.0
    }

    #[inline]
    fn exp10_hi_th() -> Self {
        39.0
    }
}
