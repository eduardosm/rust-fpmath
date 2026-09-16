use crate::double::SemiDouble;

impl crate::generic::Log10 for f32 {
    // GENERATE: log10::consts f32
    const LOG10_E_EX: SemiDouble<f32> =
        SemiDouble::with_parts(f32::from_bits(0x3EDE5000), f32::from_bits(0x38BD8A93)); // 4.342944819e-1
    const LOG10_2_HI: f32 = f32::from_bits(0x3E9A2000); // 3.010254e-1
    const LOG10_2_LO: f32 = f32::from_bits(0x369A84FC); // 4.605039e-6
}
