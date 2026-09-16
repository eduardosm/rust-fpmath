use crate::double::SemiDouble;

impl crate::generic::Reduce90Deg for f64 {
    // GENERATE: reduce_90_deg::consts f64
    const DEG_TO_RAD: f64 = f64::from_bits(0x3F91DF46A2529D39); // 1.7453292519943295e-2
    const DEG_TO_RAD_EX: SemiDouble<f64> = SemiDouble::with_parts(
        f64::from_bits(0x3F91DF46A0000000), // 1.745329238474369e-2
        f64::from_bits(0x3DE294E9C8AE0EC6), // 1.3519960527851425e-10
    );

    type SRaw = i64;
}
