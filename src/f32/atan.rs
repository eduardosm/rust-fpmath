use super::{F32Like, LikeF32};

// GENERATE: atan::consts f32
const FRAC_PI_2_HI: u32 = 0x3FC90FDA; // 1.5707963e0
const FRAC_PI_2_LO: u32 = 0x33A22169; // 7.54979e-8
const FRAC_3PI_4: u32 = 0x4016CBE4; // 2.3561945e0

impl<F: F32Like> crate::generic::Atan<LikeF32> for F {
    #[inline]
    fn frac_pi_2_hi() -> Self {
        Self::from_raw(FRAC_PI_2_HI)
    }

    #[inline]
    fn frac_pi_2_lo() -> Self {
        Self::from_raw(FRAC_PI_2_LO)
    }

    #[inline]
    fn frac_3pi_4() -> Self {
        Self::from_raw(FRAC_3PI_4)
    }

    #[inline]
    fn atan_poly(x2: Self) -> (Self, Self) {
        // GENERATE: atan::atan_poly f32 9
        const K3: u32 = 0xBEAAAA93; // -3.3333263e-1
        const K5: u32 = 0x3E4CC690; // 1.999762e-1
        const K7: u32 = 0xBE11F9D7; // -1.4255463e-1
        const K9: u32 = 0x3DDF7143; // 1.0910275e-1
        const K11: u32 = 0xBDA9C62A; // -8.28975e-2
        const K13: u32 = 0x3D65E59F; // 5.6127187e-2
        const K15: u32 = 0xBCF1397E; // -2.944636e-2
        const K17: u32 = 0x3C244237; // 1.00255525e-2
        const K19: u32 = 0xBAD2081E; // -1.6024148e-3

        let k3 = Self::from_raw(K3);
        let k5 = Self::from_raw(K5);
        let k7 = Self::from_raw(K7);
        let k9 = Self::from_raw(K9);
        let k11 = Self::from_raw(K11);
        let k13 = Self::from_raw(K13);
        let k15 = Self::from_raw(K15);
        let k17 = Self::from_raw(K17);
        let k19 = Self::from_raw(K19);

        let t = horner!(x2, x2, [k5, k7, k9, k11, k13, k15, k17, k19]);
        (k3, t)
    }
}
