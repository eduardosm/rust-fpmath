impl crate::generic::Exp10 for f32 {
    // GENERATE: exp10::consts f32
    const LOG2_10: f32 = f32::from_bits(0x40549A78); // 3.321928e0
    const LOG10_2_HI: f32 = f32::from_bits(0x3E9A2000); // 3.010254e-1
    const LOG10_2_LO: f32 = f32::from_bits(0x369A84FC); // 4.605039e-6
    const LN_10: f32 = f32::from_bits(0x40135D8E); // 2.3025851e0
    const LN_10_HI: f32 = f32::from_bits(0x40135000); // 2.3017578e0
    const LN_10_LO: f32 = f32::from_bits(0x3A58DDDB); // 8.272805e-4

    const EXP10_LO_TH: Self = -46.0;
    const EXP10_HI_TH: Self = 39.0;
}
