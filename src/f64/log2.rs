use crate::double::SemiDouble;

// GENERATE: log2::consts f64
const LOG2_E_HI: f64 = f64::from_bits(0x3FF7154760000000); // 1.4426950216293335e0
const LOG2_E_LO: f64 = f64::from_bits(0x3E54AE0BF85DDF44); // 1.9259629911266175e-8

impl crate::generic::Log2 for f64 {
    #[inline]
    fn log2_e_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(LOG2_E_HI, LOG2_E_LO)
    }
}
