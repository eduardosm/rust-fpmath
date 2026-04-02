use crate::double::SemiDouble;

// GENERATE: reduce_90_deg::consts f32
const DEG_TO_RAD: u32 = 0x3C8EFA35; // 1.7453292e-2
const DEG_TO_RAD_HI: u32 = 0x3C8EF000; // 1.7448425e-2
const DEG_TO_RAD_LO: u32 = 0x36A35129; // 4.867227e-6

impl crate::generic::Reduce90Deg for f32 {
    #[inline]
    fn deg_to_rad() -> Self {
        f32::from_bits(DEG_TO_RAD)
    }

    #[inline]
    fn deg_to_rad_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(f32::from_bits(DEG_TO_RAD_HI), f32::from_bits(DEG_TO_RAD_LO))
    }

    type SRaw = i32;
}
