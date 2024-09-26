use super::{F64Like, LikeF64};
use crate::double::SemiDouble;

// Generated with `./run-generator.sh f64::rad_to_deg::consts`
const RAD_TO_DEG: u64 = 0x404CA5DC1A63C1F8; // 5.729577951308232e1
const RAD_TO_DEG_HI: u64 = 0x404CA5DC18000000; // 5.729577922821045e1
const RAD_TO_DEG_LO: u64 = 0x3E931E0FBDC30A97; // 2.8487187165804814e-7

impl<F: F64Like> crate::generic::RadToDeg<LikeF64> for F {
    #[inline]
    fn rad_to_deg() -> Self {
        Self::from_raw(RAD_TO_DEG)
    }

    #[inline]
    fn rad_to_deg_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(Self::from_raw(RAD_TO_DEG_HI), Self::from_raw(RAD_TO_DEG_LO))
    }
}
