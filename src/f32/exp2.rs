// GENERATE: exp2::consts f32
const LN_2: f32 = f32::from_bits(0x3F317218); // 6.931472e-1

impl crate::generic::Exp2 for f32 {
    #[inline]
    fn ln_2() -> Self {
        LN_2
    }

    #[inline]
    fn exp2_lo_th() -> Self {
        -151.0
    }

    #[inline]
    fn exp2_hi_th() -> Self {
        129.0
    }
}
