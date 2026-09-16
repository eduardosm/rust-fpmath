use crate::double::SemiDouble;

impl crate::generic::Log2 for f32 {
    // GENERATE: log2::consts f32
    const LOG2_E_EX: SemiDouble<f32> = SemiDouble::with_parts(
        f32::from_bits(0x3FB8A000), // 1.4423828e0
        f32::from_bits(0x39A3B296), // 3.122284e-4
    );
}
