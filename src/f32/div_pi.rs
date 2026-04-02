use crate::double::SemiDouble;

// GENERATE: div_pi::consts f32
const FRAC_1_PI_HI: u32 = 0x3EA2F000; // 3.182373e-1
const FRAC_1_PI_LO: u32 = 0x389836E5; // 7.25815e-5

impl crate::generic::DivPi for f32 {
    #[inline]
    fn frac_1_pi_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(f32::from_bits(FRAC_1_PI_HI), f32::from_bits(FRAC_1_PI_LO))
    }
}
