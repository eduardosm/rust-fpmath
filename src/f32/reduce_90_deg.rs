use super::{F32Like, LikeF32};
use crate::double::SemiDouble;

// GENERATE: f32::reduce_90_deg::consts
const DEG_TO_RAD: u32 = 0x3C8EFA35; // 1.7453292e-2
const DEG_TO_RAD_HI: u32 = 0x3C8EF000; // 1.7448425e-2
const DEG_TO_RAD_LO: u32 = 0x36A35129; // 4.867227e-6

impl<F: F32Like> crate::generic::Reduce90Deg<LikeF32> for F {
    #[inline]
    fn deg_to_rad() -> Self {
        Self::from_raw(DEG_TO_RAD)
    }

    #[inline]
    fn deg_to_rad_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(Self::from_raw(DEG_TO_RAD_HI), Self::from_raw(DEG_TO_RAD_LO))
    }

    type SRaw = i32;
}
