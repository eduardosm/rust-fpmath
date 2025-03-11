use super::{F64Like, LikeF64};

// GENERATE: other f64::exp2::consts
const LN_2: u64 = 0x3FE62E42FEFA39EF; // 6.931471805599453e-1

impl<F: F64Like> crate::generic::Exp2<LikeF64> for F {
    #[inline]
    fn ln_2() -> Self {
        Self::from_raw(LN_2)
    }

    #[inline]
    fn exp2_lo_th() -> Self {
        Self::cast_from(-1076i16)
    }

    #[inline]
    fn exp2_hi_th() -> Self {
        Self::cast_from(1025i16)
    }
}
