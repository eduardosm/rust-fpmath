// GENERATE: exp2::consts f64
const LN_2: u64 = 0x3FE62E42FEFA39EF; // 6.931471805599453e-1

impl crate::generic::Exp2 for f64 {
    #[inline]
    fn ln_2() -> Self {
        f64::from_bits(LN_2)
    }

    #[inline]
    fn exp2_lo_th() -> Self {
        -1076.0
    }

    #[inline]
    fn exp2_hi_th() -> Self {
        1025.0
    }
}
