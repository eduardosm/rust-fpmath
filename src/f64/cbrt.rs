use crate::double::SemiDouble;

// GENERATE: cbrt::consts f64
const CBRT_2_HI: f64 = f64::from_bits(0x3FF428A2F8000000); // 1.2599210441112518e0
const CBRT_2_LO: f64 = f64::from_bits(0x3E38D728AE223DDB); // 5.783621333712523e-9
const CBRT_4_HI: f64 = f64::from_bits(0x3FF965FEA0000000); // 1.587401032447815e0
const CBRT_4_LO: f64 = f64::from_bits(0x3E54F5B8F20AC166); // 1.9520384533345454e-8

impl crate::generic::Cbrt for f64 {
    #[inline]
    fn cbrt_2_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(CBRT_2_HI, CBRT_2_LO)
    }

    #[inline]
    fn cbrt_4_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(CBRT_4_HI, CBRT_4_LO)
    }

    #[inline]
    fn exp_mod_3(e: i16) -> i8 {
        (((e + 1077) as u16) % 3) as i8
    }

    #[inline]
    fn inv_cbrt_poly(x: Self) -> Self {
        // GENERATE: cbrt::inv_cbrt_poly f64 6
        const K0: f64 = f64::from_bits(0x3FFC880B69FCA3C8); // 1.7832140102493543e0
        const K1: f64 = f64::from_bits(0xBFF92D75CD846C60); // -1.573598673631544e0
        const K2: f64 = f64::from_bits(0x3FF4116F47A08958); // 1.2542565153058458e0
        const K3: f64 = f64::from_bits(0xBFE35C39AE807700); // -6.050079735029783e-1
        const K4: f64 = f64::from_bits(0x3FC448AC781671FD); // 1.584678255429849e-1
        const K5: f64 = f64::from_bits(0xBF91C0A9E3B225C5); // -1.7336515924886848e-2

        K0 + horner!(x, x, [K1, K2, K3, K4, K5])
    }
}
