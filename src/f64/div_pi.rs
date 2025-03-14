use super::{F64Like, LikeF64};
use crate::double::SemiDouble;

// GENERATE: div_pi::consts f64
const FRAC_1_PI_HI: u64 = 0x3FD45F3068000000; // 3.1830988079309464e-1
const FRAC_1_PI_LO: u64 = 0x3E3727220A94FE14; // 5.390696036528002e-9

impl<F: F64Like> crate::generic::DivPi<LikeF64> for F {
    #[inline]
    fn frac_1_pi_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(Self::from_raw(FRAC_1_PI_HI), Self::from_raw(FRAC_1_PI_LO))
    }
}
