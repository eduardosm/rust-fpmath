use crate::double::NormDouble;

impl crate::generic::Gamma for f32 {
    #[inline]
    fn lo_th() -> Self {
        -1000.0
    }

    #[inline]
    fn hi_th() -> Self {
        1000.0
    }

    #[inline]
    fn th_1() -> Self {
        1.2
    }

    #[inline]
    fn th_2() -> Self {
        2.3
    }

    #[inline]
    fn th_3() -> Self {
        7.0
    }

    const POLY_OFF: u8 = 3;

    #[inline]
    fn half_ln_2_pi() -> NormDouble<Self> {
        // GENERATE: gamma::consts f32
        const HALF_LN_2_PI_HI: f32 = f32::from_bits(0x3F6B3F8E); // 9.189385e-1
        const HALF_LN_2_PI_LO: f32 = f32::from_bits(0x32864BEB); // 1.5634177e-8

        NormDouble::with_parts(HALF_LN_2_PI_HI, HALF_LN_2_PI_LO)
    }

    #[inline]
    fn lgamma_poly_1(x: Self) -> (Self, Self, Self, Self) {
        // GENERATE: gamma::lgamma_poly f32 12 1 0.5 1.2001
        const K1: f32 = f32::from_bits(0xBF13C468); // -5.772157e-1
        const K2: f32 = f32::from_bits(0x3F528D34); // 8.224671e-1
        const K3: f32 = f32::from_bits(0xBECD26B9); // -4.0068606e-1
        const K4: f32 = f32::from_bits(0x3E8A884D); // 2.7057114e-1
        const K5: f32 = f32::from_bits(0xBE54530B); // -2.0734803e-1
        const K6: f32 = f32::from_bits(0x3E2E5770); // 1.7025542e-1
        const K7: f32 = f32::from_bits(0xBE13CD6B); // -1.443383e-1
        const K8: f32 = f32::from_bits(0x3DD8DEBF); // 1.05893604e-1
        const K9: f32 = f32::from_bits(0xBE16B9B5); // -1.4719279e-1
        const K10: f32 = f32::from_bits(0x3E867A2E); // 2.626509e-1
        const K11: f32 = f32::from_bits(0x3F0F88B0); // 5.6067944e-1
        const K12: f32 = f32::from_bits(0x3F49DECC); // 7.8855586e-1

        let r = horner!(x, x, [K4, K5, K6, K7, K8, K9, K10, K11, K12]);
        (r, K1, K2, K3)
    }

    #[inline]
    fn lgamma_poly_2(x: Self) -> (Self, Self, Self, Self) {
        // GENERATE: gamma::lgamma_poly f32 12 2 1.1999 2.3001
        const K1: f32 = f32::from_bits(0x3ED87730); // 4.2278433e-1
        const K2: f32 = f32::from_bits(0x3EA51A66); // 3.2246703e-1
        const K3: f32 = f32::from_bits(0xBD89F002); // -6.735231e-2
        const K4: f32 = f32::from_bits(0x3CA89911); // 2.05808e-2
        const K5: f32 = f32::from_bits(0xBBF1FE68); // -7.385064e-3
        const K6: f32 = f32::from_bits(0x3B3D7FFF); // 2.8915403e-3
        const K7: f32 = f32::from_bits(0xBA9D2771); // -1.1989904e-3
        const K8: f32 = f32::from_bits(0x39FF98F1); // 4.875134e-4
        const K9: f32 = f32::from_bits(0xB9632A1C); // -2.1664094e-4
        const K10: f32 = f32::from_bits(0x3967F98E); // 2.2122843e-4
        const K11: f32 = f32::from_bits(0x3920E4F8); // 1.5344087e-4
        const K12: f32 = f32::from_bits(0x390EFD91); // 1.3636636e-4

        let r = horner!(x, x, [K4, K5, K6, K7, K8, K9, K10, K11, K12]);
        (r, K1, K2, K3)
    }

    #[inline]
    fn special_poly(x: Self) -> Self {
        // GENERATE: gamma::special_poly f32 8 2.3
        const K0: f32 = f32::from_bits(0x3DAAAAAB); // 8.3333336e-2
        const K1: f32 = f32::from_bits(0x3B638E3A); // 3.4722225e-3
        const K2: f32 = f32::from_bits(0xBB2FB999); // -2.6813506e-3
        const K3: f32 = f32::from_bits(0xB96FF7BD); // -2.2885106e-4
        const K4: f32 = f32::from_bits(0x3A4B568E); // 7.756733e-4
        const K5: f32 = f32::from_bits(0x390C8CD7); // 1.3403907e-4
        const K6: f32 = f32::from_bits(0xBA68C21E); // -8.8790234e-4
        const K7: f32 = f32::from_bits(0x3A43D4AB); // 7.4703497e-4
        const K8: f32 = f32::from_bits(0xB96EDAFB); // -2.2779025e-4

        K0 + horner!(x, x, [K1, K2, K3, K4, K5, K6, K7, K8])
    }
}
