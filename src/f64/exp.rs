use super::{F64Like, LikeF64};

// GENERATE: other f64::exp::consts
const LOG2_E: u64 = 0x3FF71547652B82FE; // 1.4426950408889634e0
const LN_2_HI: u64 = 0x3FE62E42F8000000; // 6.931471675634384e-1
const LN_2_LO: u64 = 0x3E4BE8E7BCD5E4F2; // 1.2996506893889889e-8

impl<F: F64Like> crate::generic::Exp<LikeF64> for F {
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
        Self::cast_from(-746i16)
    }

    #[inline]
    fn exp_hi_th() -> Self {
        Self::cast_from(710i16)
    }

    #[inline]
    fn exp_m1_lo_th() -> Self {
        Self::cast_from(-709i16)
    }

    #[inline]
    fn exp_m1_hi_th() -> Self {
        Self::cast_from(710i16)
    }

    #[inline]
    fn exp_special_poly(x2: Self) -> Self {
        // GENERATE: other f64::exp::exp_special_poly
        const K2: u64 = 0xBFC555555555553E; // -1.6666666666666602e-1
        const K4: u64 = 0x3F66C16C16BEBD5F; // 2.777777777701537e-3
        const K6: u64 = 0xBF11566AAF25D9EF; // -6.613756321436464e-5
        const K8: u64 = 0x3EBBBD41C686BD9E; // 1.6533902230489546e-6
        const K10: u64 = 0xBE663769EA757FAD; // -4.1381381261527065e-8

        let k2 = Self::from_raw(K2);
        let k4 = Self::from_raw(K4);
        let k6 = Self::from_raw(K6);
        let k8 = Self::from_raw(K8);
        let k10 = Self::from_raw(K10);

        horner!(x2, x2, [k2, k4, k6, k8, k10])
    }

    #[inline]
    fn exp_m1_special_poly(x2: Self) -> Self {
        // GENERATE: other f64::exp::exp_m1_special_poly
        const K2: u64 = 0xBF911111111110F5; // -1.666666666666657e-2
        const K4: u64 = 0x3F3A01A019FE5D4F; // 3.9682539681381175e-4
        const K6: u64 = 0xBEE4CE199EC6C92F; // -9.920634476444187e-6
        const K8: u64 = 0x3E90CFCAADBE80F9; // 2.505136487088992e-7
        const K10: u64 = 0xBE3AFDDC383D6C37; // -6.2844812720443466e-9

        let k2 = Self::from_raw(K2);
        let k4 = Self::from_raw(K4);
        let k6 = Self::from_raw(K6);
        let k8 = Self::from_raw(K8);
        let k10 = Self::from_raw(K10);

        F::one() + horner!(x2, x2, [k2, k4, k6, k8, k10])
    }
}
