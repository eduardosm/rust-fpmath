// GENERATE: exp::consts f32
const LOG2_E: f32 = f32::from_bits(0x3FB8AA3B); // 1.442695e0
const LN_2_HI: f32 = f32::from_bits(0x3F317000); // 6.9311523e-1
const LN_2_LO: f32 = f32::from_bits(0x3805FDF4); // 3.1946183e-5

impl crate::generic::Exp for f32 {
    #[inline]
    fn log2_e() -> Self {
        LOG2_E
    }

    #[inline]
    fn ln_2_hi() -> Self {
        LN_2_HI
    }

    #[inline]
    fn ln_2_lo() -> Self {
        LN_2_LO
    }

    #[inline]
    fn exp_lo_th() -> Self {
        -104.0
    }

    #[inline]
    fn exp_hi_th() -> Self {
        89.0
    }

    #[inline]
    fn exp_m1_lo_th() -> Self {
        -88.0
    }

    #[inline]
    fn exp_m1_hi_th() -> Self {
        89.0
    }

    #[inline]
    fn exp_special_poly(x2: Self) -> Self {
        // GENERATE: exp::exp_special_poly f32 2
        const K2: f32 = f32::from_bits(0xBE2AAA8F); // -1.6666625e-1
        const K4: f32 = f32::from_bits(0x3B35526E); // 2.766754e-3

        horner!(x2, x2, [K2, K4])
    }

    #[inline]
    fn exp_m1_special_poly(x2: Self) -> Self {
        // GENERATE: exp::exp_m1_special_poly f32 2
        const K2: f32 = f32::from_bits(0xBC888868); // -1.6666606e-2
        const K4: f32 = f32::from_bits(0x39CF2F13); // 3.951719e-4

        1.0 + horner!(x2, x2, [K2, K4])
    }
}
