use crate::double::SemiDouble;

// GENERATE: div_pi::consts f64
const FRAC_1_PI_HI: f64 = f64::from_bits(0x3FD45F3068000000); // 3.1830988079309464e-1
const FRAC_1_PI_LO: f64 = f64::from_bits(0x3E3727220A94FE14); // 5.390696036528002e-9

impl crate::generic::DivPi for f64 {
    #[inline]
    fn frac_1_pi_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(FRAC_1_PI_HI, FRAC_1_PI_LO)
    }
}
