use crate::double::SemiDouble;

impl crate::generic::DivPi for f32 {
    // GENERATE: div_pi::consts f32
    const FRAC_1_PI_EX: SemiDouble<f32> =
        SemiDouble::with_parts(f32::from_bits(0x3EA2F000), f32::from_bits(0x389836E5)); // 3.183098862e-1
}
