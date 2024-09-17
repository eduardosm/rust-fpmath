use super::{F64Like, LikeF64};

// Generated with `./run-generator.sh f64::log::consts`
const SQRT_2: u64 = 0x3FF6A09E667F3BCD; // 1.4142135623730951e0
const LN_2_HI: u64 = 0x3FE62E42F8000000; // 6.931471675634384e-1
const LN_2_LO: u64 = 0x3E4BE8E7BCD5E4F2; // 1.2996506893889889e-8
const FRAC_2_3_HI: u64 = 0x3FE5555555555555; // 6.666666666666666e-1
const FRAC_2_3_LO: u64 = 0x3C85555555555555; // 3.700743415417188e-17
const FRAC_4_10_HI: u64 = 0x3FD9999999999999; // 3.9999999999999997e-1
const FRAC_4_10_LO: u64 = 0x3C83333333333333; // 3.3306690738754695e-17

impl<F: F64Like> crate::generic::Log<LikeF64> for F {
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
    fn frac_2_3_hi() -> Self {
        Self::from_raw(FRAC_2_3_HI)
    }

    #[inline]
    fn frac_2_3_lo() -> Self {
        Self::from_raw(FRAC_2_3_LO)
    }

    #[inline]
    fn frac_4_10_hi() -> Self {
        Self::from_raw(FRAC_4_10_HI)
    }

    #[inline]
    fn frac_4_10_lo() -> Self {
        Self::from_raw(FRAC_4_10_LO)
    }

    #[inline]
    fn log_special_poly(x: Self) -> Self {
        // Generated with `./run-generator.sh f64::log::log_special_poly`
        const K2: u64 = 0x3FE5555555555592; // 6.666666666666734e-1
        const K4: u64 = 0x3FD999999997FCEC; // 3.999999999941355e-1
        const K6: u64 = 0x3FD24924941FD4B4; // 2.8571428742663874e-1
        const K8: u64 = 0x3FCC71C51FED1022; // 2.2222198542488153e-1
        const K10: u64 = 0x3FC7466413417D8C; // 1.818356603643959e-1
        const K12: u64 = 0x3FC39A17C3B61C8B; // 1.5314003998011957e-1
        const K14: u64 = 0x3FC2F07F8A21D5D4; // 1.479644226525766e-1

        let k2 = Self::from_raw(K2);
        let k4 = Self::from_raw(K4);
        let k6 = Self::from_raw(K6);
        let k8 = Self::from_raw(K8);
        let k10 = Self::from_raw(K10);
        let k12 = Self::from_raw(K12);
        let k14 = Self::from_raw(K14);

        let x2 = x * x;
        horner!(x2, x2, [k2, k4, k6, k8, k10, k12, k14])
    }

    #[inline]
    fn log_special_poly_ex(x2: Self) -> Self {
        // Generated with `./run-generator.sh f64::log::log_special_poly_ex`
        const K6: u64 = 0x3FD24924924812EE; // 2.8571428571039703e-1
        const K8: u64 = 0x3FCC71C71F9F60BA; // 2.222222237021521e-1
        const K10: u64 = 0x3FC745CF9D840A43; // 1.8181796256256363e-1
        const K12: u64 = 0x3FC3B1C4894C673B; // 1.5386254028345e-1
        const K14: u64 = 0x3FC0FB8E8D66F0D4; // 1.3267690567398083e-1
        const K16: u64 = 0x3FC0C54412F0DDDB; // 1.3102007794235146e-1

        let k6 = Self::from_raw(K6);
        let k8 = Self::from_raw(K8);
        let k10 = Self::from_raw(K10);
        let k12 = Self::from_raw(K12);
        let k14 = Self::from_raw(K14);
        let k16 = Self::from_raw(K16);

        horner!(x2, x2, [k6, k8, k10, k12, k14, k16])
    }
}
