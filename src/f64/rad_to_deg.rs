use crate::double::SemiDouble;

impl crate::generic::RadToDeg for f64 {
    // GENERATE: rad_to_deg::consts f64
    const RAD_TO_DEG: f64 = f64::from_bits(0x404CA5DC1A63C1F8); // 5.729577951308232e1
    const RAD_TO_DEG_EX: SemiDouble<f64> = SemiDouble::with_parts(
        f64::from_bits(0x404CA5DC18000000), // 5.729577922821045e1
        f64::from_bits(0x3E931E0FBDC30A97), // 2.8487187165804814e-7
    );
}
