use super::{F64Like, LikeF64};

// Generated with `./run-generator.sh f64::div_pi::consts`
const FRAC_1_PI_HI: u64 = 0x3FD45F3068000000; // 3.1830988079309464e-1
const FRAC_1_PI_LO: u64 = 0x3E3727220A94FE14; // 5.390696036528002e-9

impl<F: F64Like> crate::generic::DivPi<LikeF64> for F {
    #[inline]
    fn frac_1_pi_hi() -> Self {
        Self::from_raw(FRAC_1_PI_HI)
    }

    #[inline]
    fn frac_1_pi_lo() -> Self {
        Self::from_raw(FRAC_1_PI_LO)
    }
}
