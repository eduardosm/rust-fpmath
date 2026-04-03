use crate::double::NormDouble;

// GENERATE: log::consts f32
const SQRT_2: f32 = f32::from_bits(0x3FB504F3); // 1.4142135e0
const LN_2_HI: f32 = f32::from_bits(0x3F317000); // 6.9311523e-1
const LN_2_LO: f32 = f32::from_bits(0x3805FDF4); // 3.1946183e-5
const FRAC_2_3_HI: f32 = f32::from_bits(0x3F2AAAAA); // 6.666666e-1
const FRAC_2_3_LO: f32 = f32::from_bits(0x332AAAAB); // 3.973643e-8
const FRAC_4_10_HI: f32 = f32::from_bits(0x3ECCCCCC); // 3.9999998e-1
const FRAC_4_10_LO: f32 = f32::from_bits(0x32CCCCCD); // 2.3841858e-8

impl crate::generic::Log for f32 {
    #[inline]
    fn sqrt_2() -> Self {
        SQRT_2
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
    fn frac_2_3_ex() -> NormDouble<Self> {
        NormDouble::with_parts(FRAC_2_3_HI, FRAC_2_3_LO)
    }

    #[inline]
    fn frac_4_10_ex() -> NormDouble<Self> {
        NormDouble::with_parts(FRAC_4_10_HI, FRAC_4_10_LO)
    }

    #[inline]
    fn log_special_poly(x: Self) -> Self {
        // GENERATE: log::log_special_poly f32 4
        const K2: f32 = f32::from_bits(0x3F2AAAAA); // 6.666666e-1
        const K4: f32 = f32::from_bits(0x3ECCCD3D); // 4.0000334e-1
        const K6: f32 = f32::from_bits(0x3E921C64); // 2.8537285e-1
        const K8: f32 = f32::from_bits(0x3E717C5F); // 2.35826e-1

        let x2 = x * x;
        horner!(x2, x2, [K2, K4, K6, K8])
    }

    #[inline]
    fn log_special_poly_ex(x2: Self) -> Self {
        // GENERATE: log::log_special_poly_ex f32 3
        const K6: f32 = f32::from_bits(0x3E92495E); // 2.85716e-1
        const K8: f32 = f32::from_bits(0x3E634F16); // 2.2198138e-1
        const K10: f32 = f32::from_bits(0x3E454D62); // 1.92678e-1

        horner!(x2, x2, [K6, K8, K10])
    }
}
