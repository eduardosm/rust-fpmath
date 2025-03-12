use super::{F32Like, LikeF32};
use crate::double::NormDouble;

impl<F: F32Like> crate::generic::Gamma<LikeF32> for F {
    #[inline]
    fn lo_th() -> Self {
        Self::cast_from(-1000i16)
    }

    #[inline]
    fn hi_th() -> Self {
        Self::cast_from(1000i16)
    }

    #[inline]
    fn th_1() -> Self {
        Self::from_raw(0x3F99999A) // 1.2
    }

    #[inline]
    fn th_2() -> Self {
        Self::from_raw(0x40133333) // 2.3
    }

    #[inline]
    fn th_3() -> Self {
        Self::from_raw(0x40E00000) // 7
    }

    const POLY_OFF: u8 = 3;

    #[inline]
    fn half_ln_2_pi() -> NormDouble<Self> {
        // GENERATE: f32::gamma::consts
        const HALF_LN_2_PI_HI: u32 = 0x3F6B3F8E; // 9.189385e-1
        const HALF_LN_2_PI_LO: u32 = 0x32864BEB; // 1.5634177e-8

        NormDouble::with_parts(
            Self::from_raw(HALF_LN_2_PI_HI),
            Self::from_raw(HALF_LN_2_PI_LO),
        )
    }

    #[inline]
    fn lgamma_poly_1(x: Self) -> (Self, Self, Self, Self) {
        // GENERATE: f32::gamma::lgamma_poly_1
        const K1: u32 = 0xBF13C468; // -5.772157e-1
        const K2: u32 = 0x3F528D34; // 8.224671e-1
        const K3: u32 = 0xBECD2724; // -4.0068924e-1
        const K4: u32 = 0x3E8A8775; // 2.705647e-1
        const K5: u32 = 0xBE541CDC; // -2.0714134e-1
        const K6: u32 = 0x3E2EF709; // 1.7086424e-1
        const K7: u32 = 0xBE18A740; // -1.4907551e-1
        const K8: u32 = 0x3DAD637B; // 8.46624e-2
        const K9: u32 = 0xBE052FA3; // -1.3006453e-1
        const K10: u32 = 0x3EFA602D; // 4.89015e-1
        const K11: u32 = 0x3F7C4B16; // 9.855207e-1
        const K12: u32 = 0x3F85D08F; // 1.0454272e0

        let k1 = Self::from_raw(K1);
        let k2 = Self::from_raw(K2);
        let k3 = Self::from_raw(K3);
        let k4 = Self::from_raw(K4);
        let k5 = Self::from_raw(K5);
        let k6 = Self::from_raw(K6);
        let k7 = Self::from_raw(K7);
        let k8 = Self::from_raw(K8);
        let k9 = Self::from_raw(K9);
        let k10 = Self::from_raw(K10);
        let k11 = Self::from_raw(K11);
        let k12 = Self::from_raw(K12);

        let r = horner!(x, x, [k4, k5, k6, k7, k8, k9, k10, k11, k12]);
        (r, k1, k2, k3)
    }

    #[inline]
    fn lgamma_poly_2(x: Self) -> (Self, Self, Self, Self) {
        // GENERATE: f32::gamma::lgamma_poly_2
        const K1: u32 = 0x3ED87730; // 4.2278433e-1
        const K2: u32 = 0x3EA51A66; // 3.2246703e-1
        const K3: u32 = 0xBD89F004; // -6.7352325e-2
        const K4: u32 = 0x3CA89909; // 2.0580785e-2
        const K5: u32 = 0xBBF1FC4B; // -7.384812e-3
        const K6: u32 = 0x3B3D8928; // 2.8920863e-3
        const K7: u32 = 0xBA9D5730; // -1.2004133e-3
        const K8: u32 = 0x39FCDFA6; // 4.8231817e-4
        const K9: u32 = 0xB96522B2; // -2.1852067e-4
        const K10: u32 = 0x3971F78E; // 2.3075772e-4
        const K11: u32 = 0x392E52A9; // 1.6624726e-4
        const K12: u32 = 0x391420D1; // 1.4126605e-4

        let k1 = Self::from_raw(K1);
        let k2 = Self::from_raw(K2);
        let k3 = Self::from_raw(K3);
        let k4 = Self::from_raw(K4);
        let k5 = Self::from_raw(K5);
        let k6 = Self::from_raw(K6);
        let k7 = Self::from_raw(K7);
        let k8 = Self::from_raw(K8);
        let k9 = Self::from_raw(K9);
        let k10 = Self::from_raw(K10);
        let k11 = Self::from_raw(K11);
        let k12 = Self::from_raw(K12);

        let r = horner!(x, x, [k4, k5, k6, k7, k8, k9, k10, k11, k12]);
        (r, k1, k2, k3)
    }

    #[inline]
    fn special_poly(x: Self) -> Self {
        // GENERATE: f32::gamma::special_poly
        const K0: u32 = 0x3DAAAAAB; // 8.3333336e-2
        const K1: u32 = 0x3B638E3A; // 3.4722225e-3
        const K2: u32 = 0xBB2FB999; // -2.6813506e-3
        const K3: u32 = 0xB96FF798; // -2.2885052e-4
        const K4: u32 = 0x3A4B562E; // 7.756677e-4
        const K5: u32 = 0x390C9519; // 1.3406984e-4
        const K6: u32 = 0xBA68C85B; // -8.879953e-4
        const K7: u32 = 0x3A43DE66; // 7.4717996e-4
        const K8: u32 = 0xB96EF373; // -2.278814e-4

        let k0 = Self::from_raw(K0);
        let k1 = Self::from_raw(K1);
        let k2 = Self::from_raw(K2);
        let k3 = Self::from_raw(K3);
        let k4 = Self::from_raw(K4);
        let k5 = Self::from_raw(K5);
        let k6 = Self::from_raw(K6);
        let k7 = Self::from_raw(K7);
        let k8 = Self::from_raw(K8);

        k0 + horner!(x, x, [k1, k2, k3, k4, k5, k6, k7, k8])
    }
}
