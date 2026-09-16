use crate::double::SemiDouble;

impl crate::generic::ReduceHalfMulPi for f64 {
    // GENERATE: reduce_half_mul_pi::consts f64
    const PI_EX: SemiDouble<f64> = SemiDouble::with_parts(
        f64::from_bits(0x400921FB50000000), // 3.1415926218032837e0
        f64::from_bits(0x3E6110B4611A6263), // 3.178650954705639e-8
    );
}
