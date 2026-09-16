impl crate::generic::Exp2 for f64 {
    // GENERATE: exp2::consts f64
    const LN_2: f64 = f64::from_bits(0x3FE62E42FEFA39EF); // 6.931471805599453e-1

    const EXP2_LO_TH: Self = -1076.0;
    const EXP2_HI_TH: Self = 1025.0;
}
