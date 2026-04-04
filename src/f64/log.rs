use crate::double::NormDouble;

// GENERATE: ln::consts f64
const SQRT_2: f64 = f64::from_bits(0x3FF6A09E667F3BCD); // 1.4142135623730951e0
const LN_2_HI: f64 = f64::from_bits(0x3FE62E42F8000000); // 6.931471675634384e-1
const LN_2_LO: f64 = f64::from_bits(0x3E4BE8E7BCD5E4F2); // 1.2996506893889889e-8
const FRAC_2_3_HI: f64 = f64::from_bits(0x3FE5555555555555); // 6.666666666666666e-1
const FRAC_2_3_LO: f64 = f64::from_bits(0x3C85555555555555); // 3.700743415417188e-17
const FRAC_4_10_HI: f64 = f64::from_bits(0x3FD9999999999999); // 3.9999999999999997e-1
const FRAC_4_10_LO: f64 = f64::from_bits(0x3C83333333333333); // 3.3306690738754695e-17

impl crate::generic::Ln for f64 {
    #[inline]
    fn sqrt_2() -> Self {
        SQRT_2
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
    fn frac_2_3_ex() -> NormDouble<Self> {
        NormDouble::with_parts(FRAC_2_3_HI, FRAC_2_3_LO)
    }

    #[inline]
    fn frac_4_10_ex() -> NormDouble<Self> {
        NormDouble::with_parts(FRAC_4_10_HI, FRAC_4_10_LO)
    }

    #[inline]
    fn ln_special_poly(x: Self) -> Self {
        // GENERATE: ln::ln_special_poly f64 7
        const K2: f64 = f64::from_bits(0x3FE5555555555592); // 6.666666666666734e-1
        const K4: f64 = f64::from_bits(0x3FD999999997FCEC); // 3.999999999941355e-1
        const K6: f64 = f64::from_bits(0x3FD24924941FD4A8); // 2.8571428742663807e-1
        const K8: f64 = f64::from_bits(0x3FCC71C51FED1917); // 2.2222198542494517e-1
        const K10: f64 = f64::from_bits(0x3FC74664133FF5ED); // 1.8183566036161328e-1
        const K12: f64 = f64::from_bits(0x3FC39A17C3D56A7B); // 1.5314004003706203e-1
        const K14: f64 = f64::from_bits(0x3FC2F07F89345EF0); // 1.4796442222063133e-1

        let x2 = x * x;
        horner!(x2, x2, [K2, K4, K6, K8, K10, K12, K14])
    }

    #[inline]
    fn ln_special_poly_ex(x2: Self) -> Self {
        // GENERATE: ln::ln_special_poly_ex f64 6
        const K6: f64 = f64::from_bits(0x3FD24924924812EE); // 2.8571428571039703e-1
        const K8: f64 = f64::from_bits(0x3FCC71C71F9F60BA); // 2.222222237021521e-1
        const K10: f64 = f64::from_bits(0x3FC745CF9D840A43); // 1.8181796256256363e-1
        const K12: f64 = f64::from_bits(0x3FC3B1C4894C673B); // 1.5386254028345e-1
        const K14: f64 = f64::from_bits(0x3FC0FB8E8D66F0D4); // 1.3267690567398083e-1
        const K16: f64 = f64::from_bits(0x3FC0C54412F0DDDB); // 1.3102007794235146e-1

        horner!(x2, x2, [K6, K8, K10, K12, K14, K16])
    }
}
