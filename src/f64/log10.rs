use crate::double::SemiDouble;

// GENERATE: log10::consts f64
const LOG10_E_HI: f64 = f64::from_bits(0x3FDBCB7B10000000); // 4.342944771051407e-1
const LOG10_E_LO: f64 = f64::from_bits(0x3E349B9438CA9AAE); // 4.798111141615973e-9
const LOG10_2_HI: f64 = f64::from_bits(0x3FD3441350000000); // 3.010299950838089e-1
const LOG10_2_LO: f64 = f64::from_bits(0x3E03EF3FDE623E25); // 5.801722962879576e-10

impl crate::generic::Log10 for f64 {
    #[inline]
    fn log10_e_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(LOG10_E_HI, LOG10_E_LO)
    }

    #[inline]
    fn log10_2_hi() -> Self {
        LOG10_2_HI
    }

    #[inline]
    fn log10_2_lo() -> Self {
        LOG10_2_LO
    }
}
