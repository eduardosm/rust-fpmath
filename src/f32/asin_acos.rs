use super::{F32Like, LikeF32};
use crate::double::NormDouble;

// GENERATE: asin_acos::consts f32
const FRAC_PI_2_HI: u32 = 0x3FC90FDA; // 1.5707963e0
const FRAC_PI_2_LO: u32 = 0x33A22169; // 7.54979e-8

impl<F: F32Like> crate::generic::AsinAcos<LikeF32> for F {
    #[inline]
    fn frac_pi_2_ex() -> NormDouble<Self> {
        NormDouble::with_parts(Self::from_raw(FRAC_PI_2_HI), Self::from_raw(FRAC_PI_2_LO))
    }

    #[inline]
    fn asin_poly(x2: Self) -> Self {
        // GENERATE: asin_acos::asin_poly f32 5
        const K0: u32 = 0x3E2AAB15; // 1.6666825e-1
        const K2: u32 = 0x3D99749A; // 7.492943e-2
        const K4: u32 = 0x3D3B48FF; // 4.572391e-2
        const K6: u32 = 0x3CBCF147; // 2.3064269e-2
        const K8: u32 = 0x3D33CAF4; // 4.3894723e-2

        let k0 = Self::from_raw(K0);
        let k2 = Self::from_raw(K2);
        let k4 = Self::from_raw(K4);
        let k6 = Self::from_raw(K6);
        let k8 = Self::from_raw(K8);

        k0 + horner!(x2, x2, [k2, k4, k6, k8])
    }
}
