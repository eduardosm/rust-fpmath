impl crate::generic::Exp10 for f64 {
    // GENERATE: exp10::consts f64
    const LOG2_10: f64 = f64::from_bits(0x400A934F0979A371); // 3.321928094887362e0
    const LOG10_2_HI: f64 = f64::from_bits(0x3FD3441350000000); // 3.010299950838089e-1
    const LOG10_2_LO: f64 = f64::from_bits(0x3E03EF3FDE623E25); // 5.801722962879576e-10
    const LN_10: f64 = f64::from_bits(0x40026BB1BBB55516); // 2.302585092994046e0
    const LN_10_HI: f64 = f64::from_bits(0x40026BB1B8000000); // 2.3025850653648376e0
    const LN_10_LO: f64 = f64::from_bits(0x3E5DAAA8AC16EA57); // 2.7629208037533617e-8

    const EXP10_LO_TH: Self = -324.0;
    const EXP10_HI_TH: Self = 309.0;
}
