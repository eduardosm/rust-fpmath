use crate::double::SemiDouble;

impl crate::generic::ReduceHalfMulPi for f32 {
    // GENERATE: reduce_half_mul_pi::consts f32
    const PI_EX: SemiDouble<f32> =
        SemiDouble::with_parts(f32::from_bits(0x40490000), f32::from_bits(0x3A7DAA22)); // 3.141592654e0
}
