use crate::double::SemiDouble;

impl crate::generic::Reduce90Deg for f32 {
    // GENERATE: reduce_90_deg::consts f32
    const DEG_TO_RAD: f32 = f32::from_bits(0x3C8EFA35); // 1.7453292e-2
    const DEG_TO_RAD_EX: SemiDouble<f32> =
        SemiDouble::with_parts(f32::from_bits(0x3C8EF000), f32::from_bits(0x36A35129)); // 1.745329252e-2

    type SRaw = i32;
}
