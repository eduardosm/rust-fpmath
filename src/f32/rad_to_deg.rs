use crate::double::SemiDouble;

// GENERATE: rad_to_deg::consts f32
const RAD_TO_DEG: f32 = f32::from_bits(0x42652EE1); // 5.729578e1
const RAD_TO_DEG_HI: f32 = f32::from_bits(0x42652000); // 5.728125e1
const RAD_TO_DEG_LO: f32 = f32::from_bits(0x3C6E0D32); // 1.4529513e-2

impl crate::generic::RadToDeg for f32 {
    #[inline]
    fn rad_to_deg() -> Self {
        RAD_TO_DEG
    }

    #[inline]
    fn rad_to_deg_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(RAD_TO_DEG_HI, RAD_TO_DEG_LO)
    }
}
