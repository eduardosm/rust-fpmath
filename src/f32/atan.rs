// GENERATE: atan::consts f32
const FRAC_PI_2_HI: f32 = f32::from_bits(0x3FC90FDA); // 1.5707963e0
const FRAC_PI_2_LO: f32 = f32::from_bits(0x33A22169); // 7.54979e-8
const FRAC_3PI_4: f32 = f32::from_bits(0x4016CBE4); // 2.3561945e0

impl crate::generic::Atan for f32 {
    #[inline]
    fn frac_pi_2_hi() -> Self {
        FRAC_PI_2_HI
    }

    #[inline]
    fn frac_pi_2_lo() -> Self {
        FRAC_PI_2_LO
    }

    #[inline]
    fn frac_3pi_4() -> Self {
        FRAC_3PI_4
    }

    #[inline]
    fn atan_poly(x2: Self) -> (Self, Self) {
        // GENERATE: atan::atan_poly f32 9
        const K3: f32 = f32::from_bits(0xBEAAAA93); // -3.3333263e-1
        const K5: f32 = f32::from_bits(0x3E4CC690); // 1.999762e-1
        const K7: f32 = f32::from_bits(0xBE11F9D7); // -1.4255463e-1
        const K9: f32 = f32::from_bits(0x3DDF7143); // 1.0910275e-1
        const K11: f32 = f32::from_bits(0xBDA9C62A); // -8.28975e-2
        const K13: f32 = f32::from_bits(0x3D65E59F); // 5.6127187e-2
        const K15: f32 = f32::from_bits(0xBCF1397E); // -2.944636e-2
        const K17: f32 = f32::from_bits(0x3C244237); // 1.00255525e-2
        const K19: f32 = f32::from_bits(0xBAD2081E); // -1.6024148e-3

        let t = horner!(x2, x2, [K5, K7, K9, K11, K13, K15, K17, K19]);
        (K3, t)
    }
}
