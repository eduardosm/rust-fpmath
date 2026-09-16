use crate::double::SemiDouble;

impl crate::generic::Log2 for f64 {
    // GENERATE: log2::consts f64
    const LOG2_E_EX: SemiDouble<f64> = SemiDouble::with_parts(
        f64::from_bits(0x3FF7154760000000),
        f64::from_bits(0x3E54AE0BF85DDF44),
    ); // 1.4426950408889634073599e0
}
