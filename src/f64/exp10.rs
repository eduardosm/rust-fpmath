use super::{F64Like, LikeF64};

// GENERATE: f64::exp10::consts
const LOG2_10: u64 = 0x400A934F0979A371; // 3.321928094887362e0
const LOG10_2_HI: u64 = 0x3FD3441350000000; // 3.010299950838089e-1
const LOG10_2_LO: u64 = 0x3E03EF3FDE623E25; // 5.801722962879576e-10
const LN_10: u64 = 0x40026BB1BBB55516; // 2.302585092994046e0
const LN_10_HI: u64 = 0x40026BB1B8000000; // 2.3025850653648376e0
const LN_10_LO: u64 = 0x3E5DAAA8AC16EA57; // 2.7629208037533617e-8

impl<F: F64Like> crate::generic::Exp10<LikeF64> for F {
    #[inline]
    fn log2_10() -> Self {
        Self::from_raw(LOG2_10)
    }

    #[inline]
    fn log10_2_hi() -> Self {
        Self::from_raw(LOG10_2_HI)
    }

    #[inline]
    fn log10_2_lo() -> Self {
        Self::from_raw(LOG10_2_LO)
    }

    #[inline]
    fn ln_10() -> Self {
        Self::from_raw(LN_10)
    }

    #[inline]
    fn ln_10_hi() -> Self {
        Self::from_raw(LN_10_HI)
    }

    #[inline]
    fn ln_10_lo() -> Self {
        Self::from_raw(LN_10_LO)
    }

    #[inline]
    fn exp10_lo_th() -> Self {
        Self::cast_from(-324i16)
    }

    #[inline]
    fn exp10_hi_th() -> Self {
        Self::cast_from(309i16)
    }
}
