use crate::double::SemiDouble;

// GENERATE: rad_to_deg::consts f32
const RAD_TO_DEG: u32 = 0x42652EE1; // 5.729578e1
const RAD_TO_DEG_HI: u32 = 0x42652000; // 5.728125e1
const RAD_TO_DEG_LO: u32 = 0x3C6E0D32; // 1.4529513e-2

impl crate::generic::RadToDeg for f32 {
    #[inline]
    fn rad_to_deg() -> Self {
        f32::from_bits(RAD_TO_DEG)
    }

    #[inline]
    fn rad_to_deg_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(f32::from_bits(RAD_TO_DEG_HI), f32::from_bits(RAD_TO_DEG_LO))
    }
}
