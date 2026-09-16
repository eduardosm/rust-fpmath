impl crate::generic::Exp2 for f32 {
    // GENERATE: exp2::consts f32
    const LN_2: f32 = f32::from_bits(0x3F317218); // 6.931472e-1

    const EXP2_LO_TH: Self = -151.0;
    const EXP2_HI_TH: Self = 129.0;
}
