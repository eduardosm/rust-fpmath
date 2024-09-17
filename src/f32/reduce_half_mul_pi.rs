use super::{F32Like, LikeF32};

// Generated with `./run-generator.sh f32::reduce_half_mul_pi::consts`
const PI: u32 = 0x40490FDB; // 3.1415927e0
const PI_HI: u32 = 0x40490000; // 3.140625e0
const PI_LO: u32 = 0x3A7DAA22; // 9.676536e-4

impl<F: F32Like> crate::generic::ReduceHalfMulPi<LikeF32> for F {
    #[inline]
    fn pi() -> Self {
        Self::from_raw(PI)
    }

    #[inline]
    fn pi_hi() -> Self {
        Self::from_raw(PI_HI)
    }

    #[inline]
    fn pi_lo() -> Self {
        Self::from_raw(PI_LO)
    }
}
