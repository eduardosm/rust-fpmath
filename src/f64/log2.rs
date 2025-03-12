use super::{F64Like, LikeF64};
use crate::double::SemiDouble;

// GENERATE: f64::log2::consts
const LOG2_E_HI: u64 = 0x3FF7154760000000; // 1.4426950216293335e0
const LOG2_E_LO: u64 = 0x3E54AE0BF85DDF44; // 1.9259629911266175e-8

impl<F: F64Like> crate::generic::Log2<LikeF64> for F {
    #[inline]
    fn log2_e_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(Self::from_raw(LOG2_E_HI), Self::from_raw(LOG2_E_LO))
    }
}
