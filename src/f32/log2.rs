use crate::double::SemiDouble;

// GENERATE: log2::consts f32
const LOG2_E_HI: f32 = f32::from_bits(0x3FB8A000); // 1.4423828e0
const LOG2_E_LO: f32 = f32::from_bits(0x39A3B296); // 3.122284e-4

impl crate::generic::Log2 for f32 {
    #[inline]
    fn log2_e_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(LOG2_E_HI, LOG2_E_LO)
    }
}
