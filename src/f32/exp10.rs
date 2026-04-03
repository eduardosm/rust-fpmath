// GENERATE: exp10::consts f32
const LOG2_10: f32 = f32::from_bits(0x40549A78); // 3.321928e0
const LOG10_2_HI: f32 = f32::from_bits(0x3E9A2000); // 3.010254e-1
const LOG10_2_LO: f32 = f32::from_bits(0x369A84FC); // 4.605039e-6
const LN_10: f32 = f32::from_bits(0x40135D8E); // 2.3025851e0
const LN_10_HI: f32 = f32::from_bits(0x40135000); // 2.3017578e0
const LN_10_LO: f32 = f32::from_bits(0x3A58DDDB); // 8.272805e-4

impl crate::generic::Exp10 for f32 {
    #[inline]
    fn log2_10() -> Self {
        LOG2_10
    }

    #[inline]
    fn log10_2_hi() -> Self {
        LOG10_2_HI
    }

    #[inline]
    fn log10_2_lo() -> Self {
        LOG10_2_LO
    }

    #[inline]
    fn ln_10() -> Self {
        LN_10
    }

    #[inline]
    fn ln_10_hi() -> Self {
        LN_10_HI
    }

    #[inline]
    fn ln_10_lo() -> Self {
        LN_10_LO
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
