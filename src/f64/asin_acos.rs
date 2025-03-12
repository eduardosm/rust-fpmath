use super::{F64Like, LikeF64};
use crate::double::NormDouble;

// GENERATE: f64::asin_acos::consts
const FRAC_PI_2_HI: u64 = 0x3FF921FB54442D18; // 1.5707963267948966e0
const FRAC_PI_2_LO: u64 = 0x3C91A62633145C07; // 6.123233995736766e-17

impl<F: F64Like> crate::generic::AsinAcos<LikeF64> for F {
    #[inline]
    fn frac_pi_2_ex() -> NormDouble<Self> {
        NormDouble::with_parts(Self::from_raw(FRAC_PI_2_HI), Self::from_raw(FRAC_PI_2_LO))
    }

    #[inline]
    fn asin_poly(x2: Self) -> Self {
        // GENERATE: f64::asin_acos::asin_poly
        const K0: u64 = 0x3FC55555555555D2; // 1.6666666666667013e-1
        const K2: u64 = 0x3FB3333333324C2E; // 7.499999999917925e-2
        const K4: u64 = 0x3FA6DB6DB77D26B9; // 4.464285721640011e-2
        const K6: u64 = 0x3F9F1C718B74D800; // 3.0381940972084465e-2
        const K8: u64 = 0x3F96E8C0DAD01AA8; // 2.2372258525161698e-2
        const K10: u64 = 0x3F91C46F55A2C90D; // 1.73509021776193e-2
        const K12: u64 = 0x3F8CA61564004848; // 1.3988654245654555e-2
        const K14: u64 = 0x3F873909EF7CF228; // 1.133926164731587e-2
        const K16: u64 = 0x3F86B8816A8DD57F; // 1.1094103874465187e-2
        const K18: u64 = 0x3F652E98FB7458AB; // 2.5856960234121504e-3
        const K20: u64 = 0x3F98E2317EB55DFA; // 2.4300359114333127e-2
        const K22: u64 = 0xBF9917F6F770485D; // -2.450548062558543e-2
        const K24: u64 = 0x3FA1C880C9A9ACF3; // 3.4732842080154834e-2

        let k0 = Self::from_raw(K0);
        let k2 = Self::from_raw(K2);
        let k4 = Self::from_raw(K4);
        let k6 = Self::from_raw(K6);
        let k8 = Self::from_raw(K8);
        let k10 = Self::from_raw(K10);
        let k12 = Self::from_raw(K12);
        let k14 = Self::from_raw(K14);
        let k16 = Self::from_raw(K16);
        let k18 = Self::from_raw(K18);
        let k20 = Self::from_raw(K20);
        let k22 = Self::from_raw(K22);
        let k24 = Self::from_raw(K24);

        k0 + horner!(
            x2,
            x2,
            [k2, k4, k6, k8, k10, k12, k14, k16, k18, k20, k22, k24]
        )
    }
}
