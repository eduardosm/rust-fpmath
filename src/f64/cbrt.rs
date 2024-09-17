use super::{F64Like, LikeF64};

// Generated with `./run-generator.sh f64::cbrt::consts`
const CBRT_2_HI: u64 = 0x3FF428A2F8000000; // 1.2599210441112518e0
const CBRT_2_LO: u64 = 0x3E38D728AE223DDB; // 5.783621333712523e-9
const CBRT_4_HI: u64 = 0x3FF965FEA0000000; // 1.587401032447815e0
const CBRT_4_LO: u64 = 0x3E54F5B8F20AC166; // 1.9520384533345454e-8

impl<F: F64Like> crate::generic::Cbrt<LikeF64> for F {
    #[inline]
    fn cbrt_2_hi() -> Self {
        Self::from_raw(CBRT_2_HI)
    }

    #[inline]
    fn cbrt_2_lo() -> Self {
        Self::from_raw(CBRT_2_LO)
    }

    #[inline]
    fn cbrt_4_hi() -> F {
        Self::from_raw(CBRT_4_HI)
    }

    #[inline]
    fn cbrt_4_lo() -> F {
        Self::from_raw(CBRT_4_LO)
    }

    #[inline]
    fn exp_mod_3(e: i16) -> i8 {
        (((e + 1077) as u16) % 3) as i8
    }

    #[inline]
    fn inv_cbrt_poly(x: Self) -> Self {
        // Generated with `./run-generator.sh f64::cbrt::inv_cbrt_poly`
        const K0: u64 = 0x3FFC880B69FCA3C8; // 1.7832140102493543e0
        const K1: u64 = 0xBFF92D75CD846C60; // -1.573598673631544e0
        const K2: u64 = 0x3FF4116F47A08958; // 1.2542565153058458e0
        const K3: u64 = 0xBFE35C39AE807700; // -6.050079735029783e-1
        const K4: u64 = 0x3FC448AC781671FD; // 1.584678255429849e-1
        const K5: u64 = 0xBF91C0A9E3B225C5; // -1.7336515924886848e-2

        let k0 = Self::from_raw(K0);
        let k1 = Self::from_raw(K1);
        let k2 = Self::from_raw(K2);
        let k3 = Self::from_raw(K3);
        let k4 = Self::from_raw(K4);
        let k5 = Self::from_raw(K5);

        k0 + horner!(x, x, [k1, k2, k3, k4, k5])
    }
}
