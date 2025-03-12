use super::{F32Like, LikeF32};

// GENERATE: f32::exp::consts
const LOG2_E: u32 = 0x3FB8AA3B; // 1.442695e0
const LN_2_HI: u32 = 0x3F317000; // 6.9311523e-1
const LN_2_LO: u32 = 0x3805FDF4; // 3.1946183e-5

impl<F: F32Like> crate::generic::Exp<LikeF32> for F {
    #[inline]
    fn log2_e() -> Self {
        Self::from_raw(LOG2_E)
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
    fn exp_lo_th() -> Self {
        Self::cast_from(-104i16)
    }

    #[inline]
    fn exp_hi_th() -> Self {
        Self::cast_from(89i16)
    }

    #[inline]
    fn exp_m1_lo_th() -> Self {
        Self::cast_from(-88i16)
    }

    #[inline]
    fn exp_m1_hi_th() -> Self {
        Self::cast_from(89i16)
    }

    #[inline]
    fn exp_special_poly(x2: Self) -> Self {
        // GENERATE: f32::exp::exp_special_poly
        const K2: u32 = 0xBE2AAA8F; // -1.6666625e-1
        const K4: u32 = 0x3B35526E; // 2.766754e-3

        let k2 = Self::from_raw(K2);
        let k4 = Self::from_raw(K4);

        horner!(x2, x2, [k2, k4])
    }

    #[inline]
    fn exp_m1_special_poly(x2: Self) -> Self {
        // GENERATE: f32::exp::exp_m1_special_poly
        const K2: u32 = 0xBC888868; // -1.6666606e-2
        const K4: u32 = 0x39CF2F13; // 3.951719e-4

        let k2 = Self::from_raw(K2);
        let k4 = Self::from_raw(K4);

        F::one() + horner!(x2, x2, [k2, k4])
    }
}
