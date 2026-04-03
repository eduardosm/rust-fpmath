use crate::double::SemiDouble;

// GENERATE: reduce_half_mul_pi::consts f64
const PI_HI: f64 = f64::from_bits(0x400921FB50000000); // 3.1415926218032837e0
const PI_LO: f64 = f64::from_bits(0x3E6110B4611A6263); // 3.178650954705639e-8

impl crate::generic::ReduceHalfMulPi for f64 {
    #[inline]
    fn pi_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(PI_HI, PI_LO)
    }
}
