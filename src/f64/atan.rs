// GENERATE: atan::consts f64
const FRAC_PI_2_HI: f64 = f64::from_bits(0x3FF921FB54442D18); // 1.5707963267948966e0
const FRAC_PI_2_LO: f64 = f64::from_bits(0x3C91A62633145C07); // 6.123233995736766e-17
const FRAC_3PI_4: f64 = f64::from_bits(0x4002D97C7F3321D2); // 2.356194490192345e0

impl crate::generic::Atan for f64 {
    #[inline]
    fn frac_pi_2_hi() -> Self {
        FRAC_PI_2_HI
    }

    #[inline]
    fn frac_pi_2_lo() -> Self {
        FRAC_PI_2_LO
    }

    #[inline]
    fn frac_3pi_4() -> Self {
        FRAC_3PI_4
    }

    #[inline]
    fn atan_poly(x2: Self) -> (Self, Self) {
        // GENERATE: atan::atan_poly f64 20
        const K3: f64 = f64::from_bits(0xBFD55555555554A8); // -3.333333333333237e-1
        const K5: f64 = f64::from_bits(0x3FC999999998DDDB); // 1.99999999998666e-1
        const K7: f64 = f64::from_bits(0xBFC2492492222261); // -1.428571427861867e-1
        const K9: f64 = f64::from_bits(0x3FBC71C713BABB10); // 1.1111110908198563e-1
        const K11: f64 = f64::from_bits(0xBFB745D0D8DBF8A8); // -9.090905470290755e-2
        const K13: f64 = f64::from_bits(0x3FB3B133B63E4546); // 7.692263793275797e-2
        const K15: f64 = f64::from_bits(0xBFB110D0F16E6DDA); // -6.666284460373309e-2
        const K17: f64 = f64::from_bits(0x3FAE1ADCFBD54AF4); // 5.87987000859779e-2
        const K19: f64 = f64::from_bits(0xBFAAE24E0347F882); // -5.2507818127183206e-2
        const K21: f64 = f64::from_bits(0x3FA822298173BE98); // 4.713563637949142e-2
        const K23: f64 = f64::from_bits(0xBFA57D989E95BD38); // -4.1973847731123504e-2
        const K25: f64 = f64::from_bits(0x3FA28B28F449C44D); // 3.621795637162304e-2
        const K27: f64 = f64::from_bits(0xBF9DF6081CE612E8); // -2.9258848918760544e-2
        const K29: f64 = f64::from_bits(0x3F95BAFB6FDFA68F); // 2.122109289476853e-2
        const K31: f64 = f64::from_bits(0xBF8B0B07B02A8E04); // -1.3204631866137191e-2
        const K33: f64 = f64::from_bits(0x3F7B84A0C8C39E9E); // 6.718280852716111e-3
        const K35: f64 = f64::from_bits(0xBF65ABAC56B5B133); // -2.6453367224147594e-3
        const K37: f64 = f64::from_bits(0x3F48893B30409082); // 7.487811094546578e-4
        const K39: f64 = f64::from_bits(0xBF21A73A6FE96DF7); // -1.346834978502712e-4
        const K41: f64 = f64::from_bits(0x3EE823E8FE5305A0); // 1.151097962759606e-5

        let t = horner!(
            x2,
            x2,
            [
                K5, K7, K9, K11, K13, K15, K17, K19, K21, K23, K25, K27, K29, K31, K33, K35, K37,
                K39, K41
            ]
        );
        (K3, t)
    }
}
