use crate::double::SemiDouble;

// GENERATE: div_pi::consts f32
const FRAC_1_PI_HI: f32 = f32::from_bits(0x3EA2F000); // 3.182373e-1
const FRAC_1_PI_LO: f32 = f32::from_bits(0x389836E5); // 7.25815e-5

impl crate::generic::DivPi for f32 {
    #[inline]
    fn frac_1_pi_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(FRAC_1_PI_HI, FRAC_1_PI_LO)
    }
}
