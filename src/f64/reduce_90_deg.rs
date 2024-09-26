use super::{F64Like, LikeF64};
use crate::double::SemiDouble;

// Generated with `./run-generator.sh f64::reduce_90_deg::consts`
const DEG_TO_RAD: u64 = 0x3F91DF46A2529D39; // 1.7453292519943295e-2
const DEG_TO_RAD_HI: u64 = 0x3F91DF46A0000000; // 1.745329238474369e-2
const DEG_TO_RAD_LO: u64 = 0x3DE294E9C8AE0EC6; // 1.3519960527851425e-10

impl<F: F64Like> crate::generic::Reduce90Deg<LikeF64> for F {
    #[inline]
    fn deg_to_rad() -> Self {
        Self::from_raw(DEG_TO_RAD)
    }

    #[inline]
    fn deg_to_rad_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(Self::from_raw(DEG_TO_RAD_HI), Self::from_raw(DEG_TO_RAD_LO))
    }

    type SRaw = i64;
}
