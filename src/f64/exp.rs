// GENERATE: exp::consts f64
const LOG2_E: f64 = f64::from_bits(0x3FF71547652B82FE); // 1.4426950408889634e0
const LN_2_HI: f64 = f64::from_bits(0x3FE62E42F8000000); // 6.931471675634384e-1
const LN_2_LO: f64 = f64::from_bits(0x3E4BE8E7BCD5E4F2); // 1.2996506893889889e-8

impl crate::generic::Exp for f64 {
    #[inline]
    fn log2_e() -> Self {
        LOG2_E
    }

    #[inline]
    fn ln_2_hi() -> Self {
        LN_2_HI
    }

    #[inline]
    fn ln_2_lo() -> Self {
        LN_2_LO
    }

    #[inline]
    fn exp_lo_th() -> Self {
        -746.0
    }

    #[inline]
    fn exp_hi_th() -> Self {
        710.0
    }

    #[inline]
    fn exp_m1_lo_th() -> Self {
        -709.0
    }

    #[inline]
    fn exp_m1_hi_th() -> Self {
        710.0
    }

    #[inline]
    fn exp_special_poly(x2: Self) -> Self {
        // GENERATE: exp::exp_special_poly f64 5
        const K2: f64 = f64::from_bits(0xBFC555555555553E); // -1.6666666666666602e-1
        const K4: f64 = f64::from_bits(0x3F66C16C16BEBD5F); // 2.777777777701537e-3
        const K6: f64 = f64::from_bits(0xBF11566AAF25D9EF); // -6.613756321436464e-5
        const K8: f64 = f64::from_bits(0x3EBBBD41C686BD9E); // 1.6533902230489546e-6
        const K10: f64 = f64::from_bits(0xBE663769EA757FAD); // -4.1381381261527065e-8

        horner!(x2, x2, [K2, K4, K6, K8, K10])
    }

    #[inline]
    fn exp_m1_special_poly(x2: Self) -> Self {
        // GENERATE: exp::exp_m1_special_poly f64 5
        const K2: f64 = f64::from_bits(0xBF911111111110F5); // -1.666666666666657e-2
        const K4: f64 = f64::from_bits(0x3F3A01A019FE5D4F); // 3.9682539681381175e-4
        const K6: f64 = f64::from_bits(0xBEE4CE199EC6C92F); // -9.920634476444187e-6
        const K8: f64 = f64::from_bits(0x3E90CFCAADBE80F9); // 2.505136487088992e-7
        const K10: f64 = f64::from_bits(0xBE3AFDDC383D6C37); // -6.2844812720443466e-9

        1.0 + horner!(x2, x2, [K2, K4, K6, K8, K10])
    }
}
