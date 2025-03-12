use super::{F64Like, LikeF64};
use crate::double::SemiDouble;

// GENERATE: f64::log10::consts
const LOG10_E_HI: u64 = 0x3FDBCB7B10000000; // 4.342944771051407e-1
const LOG10_E_LO: u64 = 0x3E349B9438CA9AAE; // 4.798111141615973e-9
const LOG10_2_HI: u64 = 0x3FD3441350000000; // 3.010299950838089e-1
const LOG10_2_LO: u64 = 0x3E03EF3FDE623E25; // 5.801722962879576e-10

impl<F: F64Like> crate::generic::Log10<LikeF64> for F {
    #[inline]
    fn log10_e_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(Self::from_raw(LOG10_E_HI), Self::from_raw(LOG10_E_LO))
    }

    #[inline]
    fn log10_2_hi() -> Self {
        Self::from_raw(LOG10_2_HI)
    }

    #[inline]
    fn log10_2_lo() -> Self {
        Self::from_raw(LOG10_2_LO)
    }
}
