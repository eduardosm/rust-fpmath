// GENERATE: atan::consts f64
const FRAC_PI_2_HI: u64 = 0x3FF921FB54442D18; // 1.5707963267948966e0
const FRAC_PI_2_LO: u64 = 0x3C91A62633145C07; // 6.123233995736766e-17
const FRAC_3PI_4: u64 = 0x4002D97C7F3321D2; // 2.356194490192345e0

impl crate::generic::Atan for f64 {
    #[inline]
    fn frac_pi_2_hi() -> Self {
        f64::from_bits(FRAC_PI_2_HI)
    }

    #[inline]
    fn frac_pi_2_lo() -> Self {
        f64::from_bits(FRAC_PI_2_LO)
    }

    #[inline]
    fn frac_3pi_4() -> Self {
        f64::from_bits(FRAC_3PI_4)
    }

    #[inline]
    fn atan_poly(x2: Self) -> (Self, Self) {
        // GENERATE: atan::atan_poly f64 20
        const K3: u64 = 0xBFD55555555554A8; // -3.333333333333237e-1
        const K5: u64 = 0x3FC999999998DDDB; // 1.99999999998666e-1
        const K7: u64 = 0xBFC2492492222261; // -1.428571427861867e-1
        const K9: u64 = 0x3FBC71C713BABB10; // 1.1111110908198563e-1
        const K11: u64 = 0xBFB745D0D8DBF8A8; // -9.090905470290755e-2
        const K13: u64 = 0x3FB3B133B63E4546; // 7.692263793275797e-2
        const K15: u64 = 0xBFB110D0F16E6DDA; // -6.666284460373309e-2
        const K17: u64 = 0x3FAE1ADCFBD54AF4; // 5.87987000859779e-2
        const K19: u64 = 0xBFAAE24E0347F882; // -5.2507818127183206e-2
        const K21: u64 = 0x3FA822298173BE98; // 4.713563637949142e-2
        const K23: u64 = 0xBFA57D989E95BD38; // -4.1973847731123504e-2
        const K25: u64 = 0x3FA28B28F449C44D; // 3.621795637162304e-2
        const K27: u64 = 0xBF9DF6081CE612E8; // -2.9258848918760544e-2
        const K29: u64 = 0x3F95BAFB6FDFA68F; // 2.122109289476853e-2
        const K31: u64 = 0xBF8B0B07B02A8E04; // -1.3204631866137191e-2
        const K33: u64 = 0x3F7B84A0C8C39E9E; // 6.718280852716111e-3
        const K35: u64 = 0xBF65ABAC56B5B133; // -2.6453367224147594e-3
        const K37: u64 = 0x3F48893B30409082; // 7.487811094546578e-4
        const K39: u64 = 0xBF21A73A6FE96DF7; // -1.346834978502712e-4
        const K41: u64 = 0x3EE823E8FE5305A0; // 1.151097962759606e-5

        let k3 = f64::from_bits(K3);
        let k5 = f64::from_bits(K5);
        let k7 = f64::from_bits(K7);
        let k9 = f64::from_bits(K9);
        let k11 = f64::from_bits(K11);
        let k13 = f64::from_bits(K13);
        let k15 = f64::from_bits(K15);
        let k17 = f64::from_bits(K17);
        let k19 = f64::from_bits(K19);
        let k21 = f64::from_bits(K21);
        let k23 = f64::from_bits(K23);
        let k25 = f64::from_bits(K25);
        let k27 = f64::from_bits(K27);
        let k29 = f64::from_bits(K29);
        let k31 = f64::from_bits(K31);
        let k33 = f64::from_bits(K33);
        let k35 = f64::from_bits(K35);
        let k37 = f64::from_bits(K37);
        let k39 = f64::from_bits(K39);
        let k41 = f64::from_bits(K41);

        let t = horner!(
            x2,
            x2,
            [
                k5, k7, k9, k11, k13, k15, k17, k19, k21, k23, k25, k27, k29, k31, k33, k35, k37,
                k39, k41
            ]
        );
        (k3, t)
    }
}
