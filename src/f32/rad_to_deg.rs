use super::{F32Like, LikeF32};
use crate::double::SemiDouble;

// GENERATE: other f32::rad_to_deg::consts
const RAD_TO_DEG: u32 = 0x42652EE1; // 5.729578e1
const RAD_TO_DEG_HI: u32 = 0x42652000; // 5.728125e1
const RAD_TO_DEG_LO: u32 = 0x3C6E0D32; // 1.4529513e-2

impl<F: F32Like> crate::generic::RadToDeg<LikeF32> for F {
    #[inline]
    fn rad_to_deg() -> Self {
        Self::from_raw(RAD_TO_DEG)
    }

    #[inline]
    fn rad_to_deg_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(Self::from_raw(RAD_TO_DEG_HI), Self::from_raw(RAD_TO_DEG_LO))
    }
}
