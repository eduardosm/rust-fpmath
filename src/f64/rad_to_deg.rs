use crate::double::SemiDouble;

// GENERATE: rad_to_deg::consts f64
const RAD_TO_DEG: u64 = 0x404CA5DC1A63C1F8; // 5.729577951308232e1
const RAD_TO_DEG_HI: u64 = 0x404CA5DC18000000; // 5.729577922821045e1
const RAD_TO_DEG_LO: u64 = 0x3E931E0FBDC30A97; // 2.8487187165804814e-7

impl crate::generic::RadToDeg for f64 {
    #[inline]
    fn rad_to_deg() -> Self {
        f64::from_bits(RAD_TO_DEG)
    }

    #[inline]
    fn rad_to_deg_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(f64::from_bits(RAD_TO_DEG_HI), f64::from_bits(RAD_TO_DEG_LO))
    }
}
