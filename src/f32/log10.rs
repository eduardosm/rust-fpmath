use crate::double::SemiDouble;

// GENERATE: log10::consts f32
const LOG10_E_HI: f32 = f32::from_bits(0x3EDE5000); // 4.342041e-1
const LOG10_E_LO: f32 = f32::from_bits(0x38BD8A93); // 9.038034e-5
const LOG10_2_HI: f32 = f32::from_bits(0x3E9A2000); // 3.010254e-1
const LOG10_2_LO: f32 = f32::from_bits(0x369A84FC); // 4.605039e-6

impl crate::generic::Log10 for f32 {
    #[inline]
    fn log10_e_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(LOG10_E_HI, LOG10_E_LO)
    }

    #[inline]
    fn log10_2_hi() -> Self {
        LOG10_2_HI
    }

    #[inline]
    fn log10_2_lo() -> Self {
        LOG10_2_LO
    }
}
