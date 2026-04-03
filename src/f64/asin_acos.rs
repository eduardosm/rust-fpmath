use crate::double::NormDouble;

// GENERATE: asin_acos::consts f64
const FRAC_PI_2_HI: f64 = f64::from_bits(0x3FF921FB54442D18); // 1.5707963267948966e0
const FRAC_PI_2_LO: f64 = f64::from_bits(0x3C91A62633145C07); // 6.123233995736766e-17

impl crate::generic::AsinAcos for f64 {
    #[inline]
    fn frac_pi_2_ex() -> NormDouble<Self> {
        NormDouble::with_parts(FRAC_PI_2_HI, FRAC_PI_2_LO)
    }

    #[inline]
    fn asin_poly(x2: Self) -> Self {
        // GENERATE: asin_acos::asin_poly f64 13
        const K0: f64 = f64::from_bits(0x3FC55555555555D2); // 1.6666666666667013e-1
        const K2: f64 = f64::from_bits(0x3FB3333333324C2E); // 7.499999999917925e-2
        const K4: f64 = f64::from_bits(0x3FA6DB6DB77D26B9); // 4.464285721640011e-2
        const K6: f64 = f64::from_bits(0x3F9F1C718B74D800); // 3.0381940972084465e-2
        const K8: f64 = f64::from_bits(0x3F96E8C0DAD01AA8); // 2.2372258525161698e-2
        const K10: f64 = f64::from_bits(0x3F91C46F55A2C90D); // 1.73509021776193e-2
        const K12: f64 = f64::from_bits(0x3F8CA61564004848); // 1.3988654245654555e-2
        const K14: f64 = f64::from_bits(0x3F873909EF7CF228); // 1.133926164731587e-2
        const K16: f64 = f64::from_bits(0x3F86B8816A8DD57F); // 1.1094103874465187e-2
        const K18: f64 = f64::from_bits(0x3F652E98FB7458AB); // 2.5856960234121504e-3
        const K20: f64 = f64::from_bits(0x3F98E2317EB55DFA); // 2.4300359114333127e-2
        const K22: f64 = f64::from_bits(0xBF9917F6F770485D); // -2.450548062558543e-2
        const K24: f64 = f64::from_bits(0x3FA1C880C9A9ACF3); // 3.4732842080154834e-2

        K0 + horner!(
            x2,
            x2,
            [K2, K4, K6, K8, K10, K12, K14, K16, K18, K20, K22, K24]
        )
    }
}
