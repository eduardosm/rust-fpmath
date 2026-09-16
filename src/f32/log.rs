use crate::double::NormDouble;

impl crate::generic::Ln for f32 {
    // GENERATE: ln::consts f32
    const SQRT_2: f32 = f32::from_bits(0x3FB504F3); // 1.4142135e0
    const LN_2_HI: f32 = f32::from_bits(0x3F317000); // 6.9311523e-1
    const LN_2_LO: f32 = f32::from_bits(0x3805FDF4); // 3.1946183e-5
    const FRAC_2_3_EX: NormDouble<f32> =
        NormDouble::with_parts(f32::from_bits(0x3F2AAAAA), f32::from_bits(0x332AAAAB)); // 6.6666666666667e-1
    const FRAC_4_10_EX: NormDouble<f32> =
        NormDouble::with_parts(f32::from_bits(0x3ECCCCCC), f32::from_bits(0x32CCCCCD)); // 4.0000000000000e-1

    #[inline]
    fn ln_special_poly(x: Self) -> Self {
        // GENERATE: ln::ln_special_poly f32 4
        const K2: f32 = f32::from_bits(0x3F2AAAAA); // 6.666666e-1
        const K4: f32 = f32::from_bits(0x3ECCCD3D); // 4.0000334e-1
        const K6: f32 = f32::from_bits(0x3E921C64); // 2.8537285e-1
        const K8: f32 = f32::from_bits(0x3E717C5F); // 2.35826e-1

        let x2 = x * x;
        horner!(x2, x2, [K2, K4, K6, K8])
    }

    #[inline]
    fn ln_special_poly_ex(x2: Self) -> Self {
        // GENERATE: ln::ln_special_poly_ex f32 3
        const K6: f32 = f32::from_bits(0x3E92495E); // 2.85716e-1
        const K8: f32 = f32::from_bits(0x3E634F16); // 2.2198138e-1
        const K10: f32 = f32::from_bits(0x3E454D62); // 1.92678e-1

        horner!(x2, x2, [K6, K8, K10])
    }
}
