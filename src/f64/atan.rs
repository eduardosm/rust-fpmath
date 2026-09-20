impl crate::generic::Atan for f64 {
    // GENERATE: atan::consts f64
    const FRAC_PI_2_HI: f64 = f64::from_bits(0x3FF921FB54442D18); // 1.5707963267948966e0
    const FRAC_PI_2_LO: f64 = f64::from_bits(0x3C91A62633145C07); // 6.123233995736766e-17
    const FRAC_3PI_4: f64 = f64::from_bits(0x4002D97C7F3321D2); // 2.356194490192345e0

    #[inline]
    fn atan_poly(x2: Self) -> (Self, Self, Self, Self) {
        // GENERATE: atan::atan_poly f64 26
        const K3: f64 = f64::from_bits(0xBFD5555555555555); // -3.333333333333333e-1
        const K5: f64 = f64::from_bits(0x3FC9999999999996); // 1.999999999999999e-1
        const K7: f64 = f64::from_bits(0xBFC2492492492381); // -1.4285714285713527e-1
        const K9: f64 = f64::from_bits(0x3FBC71C71C716302); // 1.1111111111075547e-1
        const K11: f64 = f64::from_bits(0xBFB745D1745199FC); // -9.090909089864135e-2
        const K13: f64 = f64::from_bits(0x3FB3B13B12CA37D5); // 7.692307671297212e-2
        const K15: f64 = f64::from_bits(0xBFB1111103ED3F91); // -6.666666360729169e-2
        const K17: f64 = f64::from_bits(0x3FAE1E1CFDBB58F0); // 5.882349583923652e-2
        const K19: f64 = f64::from_bits(0xBFAAF27D24A999B1); // -5.263129303493786e-2
        const K21: f64 = f64::from_bits(0x3FA8614541B1042C); // 4.761711527726961e-2
        const K23: f64 = f64::from_bits(0xBFA6416676DE8E9B); // -4.346771431963e-2
        const K25: f64 = f64::from_bits(0x3FA474B3B1B2FEEC); // 3.995286506824222e-2
        const K27: f64 = f64::from_bits(0xBFA2DFA71250E1BE); // -3.6862584100564075e-2
        const K29: f64 = f64::from_bits(0x3FA160FA58F7969F); // 3.39430078127354e-2
        const K31: f64 = f64::from_bits(0xBF9F975192DE155D); // -3.08506723299805e-2
        const K33: f64 = f64::from_bits(0x3F9BD6EC904B8537); // 2.7187057787960014e-2
        const K35: f64 = f64::from_bits(0xBF973696D2F5B81E); // -2.26691786060017e-2
        const K37: f64 = f64::from_bits(0x3F91CBE8A4EAAAC3); // 1.7379412713732344e-2
        const K39: f64 = f64::from_bits(0xBF8856F7641084CA); // -1.188462518809032e-2
        const K41: f64 = f64::from_bits(0x3F7CCAAA2879FB47); // 7.029213600275259e-3
        const K43: f64 = f64::from_bits(0xBF6C86C8857B47F6); // -3.482238429435287e-3
        const K45: f64 = f64::from_bits(0x3F56D4495CAF0D53); // 1.3933864999991677e-3
        const K47: f64 = f64::from_bits(0xBF3C2DF88464E442); // -4.2998616517625845e-4
        const K49: f64 = f64::from_bits(0x3F191016A7DAB658); // 9.560716895392599e-5
        const K51: f64 = f64::from_bits(0xBEEC7F2278351171); // -1.3588247162695591e-5
        const K53: f64 = f64::from_bits(0x3EAF09068503FCD7); // 9.24922696682232e-7

        let t = K9
            + horner!(
                x2,
                x2,
                [
                    K11, K13, K15, K17, K19, K21, K23, K25, K27, K29, K31, K33, K35, K37, K39, K41,
                    K43, K45, K47, K49, K51, K53
                ]
            );
        (K3, K5, K7, t)
    }
}
