use crate::double::SemiDouble;

impl crate::generic::RadToDeg for f32 {
    // GENERATE: rad_to_deg::consts f32
    const RAD_TO_DEG: f32 = f32::from_bits(0x42652EE1); // 5.729578e1
    const RAD_TO_DEG_EX: SemiDouble<f32> = SemiDouble::with_parts(
        f32::from_bits(0x42652000), // 5.728125e1
        f32::from_bits(0x3C6E0D32), // 1.4529513e-2
    );
}
