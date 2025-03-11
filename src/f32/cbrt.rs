use super::{F32Like, LikeF32};
use crate::double::SemiDouble;

// GENERATE: other f32::cbrt::consts
const CBRT_2_HI: u32 = 0x3FA14000; // 1.2597656e0
const CBRT_2_LO: u32 = 0x3922F98D; // 1.5542489e-4
const CBRT_4_HI: u32 = 0x3FCB2000; // 1.5869141e0
const CBRT_4_LO: u32 = 0x39FF529F; // 4.8698948e-4

impl<F: F32Like> crate::generic::Cbrt<LikeF32> for F {
    #[inline]
    fn cbrt_2_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(Self::from_raw(CBRT_2_HI), Self::from_raw(CBRT_2_LO))
    }

    #[inline]
    fn cbrt_4_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(Self::from_raw(CBRT_4_HI), Self::from_raw(CBRT_4_LO))
    }

    #[inline]
    fn exp_mod_3(e: i16) -> i8 {
        (((e + 153) as u16) % 3) as i8
    }

    #[inline]
    fn inv_cbrt_poly(x: Self) -> Self {
        // GENERATE: other f32::cbrt::inv_cbrt_poly
        const K0: u32 = 0x3FB21939; // 1.3913947e0
        const K1: u32 = 0xBEF9C752; // -4.8784882e-1
        const K2: u32 = 0x3DC257A9; // 9.489376e-2

        let k0 = Self::from_raw(K0);
        let k1 = Self::from_raw(K1);
        let k2 = Self::from_raw(K2);

        k0 + horner!(x, x, [k1, k2])
    }
}
