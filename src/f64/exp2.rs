use super::{F64Like, LikeF64};

// Generated with `./run-generator.sh f64::exp2::consts`
const LN_2: u64 = 0x3FE62E42FEFA39EF; // 6.931471805599453e-1
const LN_2_HI: u64 = 0x3FE62E42F8000000; // 6.931471675634384e-1
const LN_2_LO: u64 = 0x3E4BE8E7BCD5E4F2; // 1.2996506893889889e-8

impl<F: F64Like> crate::generic::Exp2<LikeF64> for F {
    #[inline]
    fn ln_2() -> Self {
        Self::from_raw(LN_2)
    }

    #[inline]
    fn ln_2_hi() -> Self {
        Self::from_raw(LN_2_HI)
    }

    #[inline]
    fn ln_2_lo() -> Self {
        Self::from_raw(LN_2_LO)
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
