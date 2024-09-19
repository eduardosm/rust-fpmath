use super::{F32Like, LikeF32};

// Generated with `./run-generator.sh f32::div_pi::consts`
const FRAC_1_PI_HI: u32 = 0x3EA2F000; // 3.182373e-1
const FRAC_1_PI_LO: u32 = 0x389836E5; // 7.25815e-5

impl<F: F32Like> crate::generic::DivPi<LikeF32> for F {
    #[inline]
    fn frac_1_pi_hi() -> Self {
        Self::from_raw(FRAC_1_PI_HI)
    }

    #[inline]
    fn frac_1_pi_lo() -> Self {
        Self::from_raw(FRAC_1_PI_LO)
    }
}
