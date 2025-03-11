use super::{F32Like, LikeF32};
use crate::double::NormDouble;

// GENERATE: other f32::log::consts
const SQRT_2: u32 = 0x3FB504F3; // 1.4142135e0
const LN_2_HI: u32 = 0x3F317000; // 6.9311523e-1
const LN_2_LO: u32 = 0x3805FDF4; // 3.1946183e-5
const FRAC_2_3_HI: u32 = 0x3F2AAAAA; // 6.666666e-1
const FRAC_2_3_LO: u32 = 0x332AAAAB; // 3.973643e-8
const FRAC_4_10_HI: u32 = 0x3ECCCCCC; // 3.9999998e-1
const FRAC_4_10_LO: u32 = 0x32CCCCCD; // 2.3841858e-8

impl<F: F32Like> crate::generic::Log<LikeF32> for F {
    #[inline]
    fn sqrt_2() -> Self {
        Self::from_raw(SQRT_2)
    }

    #[inline]
    fn ln_2_hi() -> Self {
        Self::from_raw(LN_2_HI)
    }

    #[inline]
    fn ln_2_lo() -> Self {
        Self::from_raw(LN_2_LO)
    }

    #[inline]
    fn frac_2_3_ex() -> NormDouble<Self> {
        NormDouble::with_parts(Self::from_raw(FRAC_2_3_HI), Self::from_raw(FRAC_2_3_LO))
    }

    #[inline]
    fn frac_4_10_ex() -> NormDouble<Self> {
        NormDouble::with_parts(Self::from_raw(FRAC_4_10_HI), Self::from_raw(FRAC_4_10_LO))
    }

    #[inline]
    fn log_special_poly(x: Self) -> Self {
        // GENERATE: other f32::log::log_special_poly
        const K2: u32 = 0x3F2AAAAA; // 6.666666e-1
        const K4: u32 = 0x3ECCCD3D; // 4.0000334e-1
        const K6: u32 = 0x3E921C64; // 2.8537285e-1
        const K8: u32 = 0x3E717C5F; // 2.35826e-1

        let k2 = Self::from_raw(K2);
        let k4 = Self::from_raw(K4);
        let k6 = Self::from_raw(K6);
        let k8 = Self::from_raw(K8);

        let x2 = x * x;
        horner!(x2, x2, [k2, k4, k6, k8])
    }

    #[inline]
    fn log_special_poly_ex(x2: Self) -> Self {
        // GENERATE: other f32::log::log_special_poly_ex
        const K6: u32 = 0x3E92495E; // 2.85716e-1
        const K8: u32 = 0x3E634F17; // 2.2198139e-1
        const K10: u32 = 0x3E454D62; // 1.92678e-1

        let k6 = Self::from_raw(K6);
        let k8 = Self::from_raw(K8);
        let k10 = Self::from_raw(K10);

        horner!(x2, x2, [k6, k8, k10])
    }
}
