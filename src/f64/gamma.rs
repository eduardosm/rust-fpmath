use super::{F64Like, LikeF64};
use crate::double::NormDouble;

impl<F: F64Like> crate::generic::Gamma<LikeF64> for F {
    #[inline]
    fn lo_th() -> Self {
        Self::cast_from(-10000i16)
    }

    #[inline]
    fn hi_th() -> Self {
        Self::cast_from(10000i16)
    }

    #[inline]
    fn th_1() -> Self {
        Self::from_raw(0x3FF199999999999A) // 1.1
    }

    #[inline]
    fn th_2() -> Self {
        Self::from_raw(0x4002666666666666) // 2.3
    }

    #[inline]
    fn th_3() -> Self {
        Self::from_raw(0x401C000000000000) // 7
    }

    const POLY_OFF: u8 = 5;

    #[inline]
    fn half_ln_2_pi() -> NormDouble<Self> {
        // GENERATE: gamma::consts f64
        const HALF_LN_2_PI_HI: u64 = 0x3FED67F1C864BEB4; // 9.189385332046727e-1
        const HALF_LN_2_PI_LO: u64 = 0x3C94D252F2400510; // 7.223936088184323e-17

        NormDouble::with_parts(
            Self::from_raw(HALF_LN_2_PI_HI),
            Self::from_raw(HALF_LN_2_PI_LO),
        )
    }

    #[inline]
    fn lgamma_poly_1(x: Self) -> (Self, Self, Self, Self) {
        // GENERATE: gamma::lgamma_poly_1 f64 26 1.1
        const K1: u64 = 0xBFE2788CFC6FB619; // -5.772156649015329e-1
        const K2: u64 = 0x3FEA51A6625307D3; // 8.224670334241132e-1
        const K3: u64 = 0xBFD9A4D55BEAB2C7; // -4.0068563438653054e-1
        const K4: u64 = 0x3FD151322AC7D624; // 2.705808084277541e-1
        const K5: u64 = 0xBFCA8B9C17AACFE0; // -2.0738555102945977e-1
        const K6: u64 = 0x3FC5B40CB105F2AD; // 1.6955717700684172e-1
        const K7: u64 = 0xBFC2703A1D222293; // -1.4404989645506258e-1
        const K8: u64 = 0x3FC010B36A674F5D; // 1.2550966926084586e-1
        const K9: u64 = 0xBFBC806804DE7F45; // -1.1133432501673919e-1
        const K10: u64 = 0x3FB9A0183CAEF329; // 1.0009910089042919e-1
        const K11: u64 = 0xBFB7487ABA16AA7F; // -9.094969790627692e-2
        const K12: u64 = 0x3FB55AC4DCD8BE68; // 8.34162749458699e-2
        const K13: u64 = 0xBFB3AAB76DDCC2E6; // -7.682367736994031e-2
        const K14: u64 = 0x3FB17E9E8AAF70B5; // 6.833830724594032e-2
        const K15: u64 = 0xBFB8148C5EB43FC0; // -9.406354248153459e-2
        const K16: u64 = 0xBF4AB531D1CBB53F; // -8.150571116779872e-4
        const K17: u64 = 0x3FDD12C0C1A8FB64; // 4.542695895396973e-1
        const K18: u64 = 0x40165DFCA8F24959; // 5.591784133708949e0
        const K19: u64 = 0x403B78CFF49C2C1E; // 2.7471923149230967e1
        const K20: u64 = 0x405658F7A41729C9; // 8.939011480581316e1
        const K21: u64 = 0x406994530189D081; // 2.046351325694741e2
        const K22: u64 = 0x407510762AE9660F; // 3.370288495175491e2
        const K23: u64 = 0x40789675FCCB8053; // 3.93403805537154e2
        const K24: u64 = 0x407377395F4EA94E; // 3.114515069077214e2
        const K25: u64 = 0x4062DBDDB10987AB; // 1.5087081195699042e2
        const K26: u64 = 0x404117D92E3A43A7; // 3.418631532521186e1

        let k1 = Self::from_raw(K1);
        let k2 = Self::from_raw(K2);
        let k3 = Self::from_raw(K3);
        let k4 = Self::from_raw(K4);
        let k5 = Self::from_raw(K5);
        let k6 = Self::from_raw(K6);
        let k7 = Self::from_raw(K7);
        let k8 = Self::from_raw(K8);
        let k9 = Self::from_raw(K9);
        let k10 = Self::from_raw(K10);
        let k11 = Self::from_raw(K11);
        let k12 = Self::from_raw(K12);
        let k13 = Self::from_raw(K13);
        let k14 = Self::from_raw(K14);
        let k15 = Self::from_raw(K15);
        let k16 = Self::from_raw(K16);
        let k17 = Self::from_raw(K17);
        let k18 = Self::from_raw(K18);
        let k19 = Self::from_raw(K19);
        let k20 = Self::from_raw(K20);
        let k21 = Self::from_raw(K21);
        let k22 = Self::from_raw(K22);
        let k23 = Self::from_raw(K23);
        let k24 = Self::from_raw(K24);
        let k25 = Self::from_raw(K25);
        let k26 = Self::from_raw(K26);

        let r = horner!(
            x,
            x,
            [
                k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18, k19, k20, k21,
                k22, k23, k24, k25, k26
            ]
        );
        (r, k1, k2, k3)
    }

    #[inline]
    fn lgamma_poly_2(x: Self) -> (Self, Self, Self, Self) {
        // GENERATE: gamma::lgamma_poly_2 f64 26 1.1 2.3
        const K1: u64 = 0x3FDB0EE6072093CE; // 4.2278433509846713e-1
        const K2: u64 = 0x3FD4A34CC4A60FA6; // 3.224670334241132e-1
        const K3: u64 = 0xBFB13E001A557607; // -6.73523010531981e-2
        const K4: u64 = 0x3F951322AC7D843A; // 2.0580808427784293e-2
        const K5: u64 = 0xBF7E404FC218F57A; // -7.385551028673882e-3
        const K6: u64 = 0x3F67ADD6EADC34D2; // 2.890510330763798e-3
        const K7: u64 = 0xBF538AC5C2BD4F12; // -1.1927539116713443e-3
        const K8: u64 = 0x3F40B36AF7DEDC08; // 5.096695237999391e-4
        const K9: u64 = 0xBF2D3FD4CD99ACA3; // -2.2315476126108955e-4
        const K10: u64 = 0x3F1A127B67287811; // 9.945753280546389e-5
        const K11: u64 = 0xBF078DE263862A3F; // -4.4926139199115155e-5
        const K12: u64 = 0x3EF580D20A43B0E4; // 2.0507054288396937e-5
        const K13: u64 = 0xBEE3CCAFCE27C2A3; // -9.441164768171395e-6
        const K14: u64 = 0x3ED257F162658760; // 4.373437639011924e-6
        const K15: u64 = 0xBEC0FEBAD845E64A; // -2.0259664685301373e-6
        const K16: u64 = 0x3EB0A769AB5E6EDC; // 9.926531396646932e-7
        const K17: u64 = 0xBE9F079F7E2B40C9; // -4.623797151106122e-7
        const K18: u64 = 0xBE5586C53F5574E9; // -2.0048066281187224e-8
        const K19: u64 = 0xBEA0F68137677AFF; // -5.055340881984814e-7
        const K20: u64 = 0x3E7B1D0666DEFDCB; // 1.0100520750257138e-7
        const K21: u64 = 0x3EB6181175C86AA4; // 1.3169060003650557e-6
        const K22: u64 = 0x3EC5B2CD6B0C00F5; // 2.5866564431488087e-6
        const K23: u64 = 0x3EC63140B7C8FE21; // 2.645539494388968e-6
        const K24: u64 = 0x3EBB94EE67D692BC; // 1.6440011728426515e-6
        const K25: u64 = 0x3EA39324890BEFAF; // 5.833737920231137e-7
        const K26: u64 = 0x3E797EF31705BA1B; // 9.497961684309852e-8

        let k1 = Self::from_raw(K1);
        let k2 = Self::from_raw(K2);
        let k3 = Self::from_raw(K3);
        let k4 = Self::from_raw(K4);
        let k5 = Self::from_raw(K5);
        let k6 = Self::from_raw(K6);
        let k7 = Self::from_raw(K7);
        let k8 = Self::from_raw(K8);
        let k9 = Self::from_raw(K9);
        let k10 = Self::from_raw(K10);
        let k11 = Self::from_raw(K11);
        let k12 = Self::from_raw(K12);
        let k13 = Self::from_raw(K13);
        let k14 = Self::from_raw(K14);
        let k15 = Self::from_raw(K15);
        let k16 = Self::from_raw(K16);
        let k17 = Self::from_raw(K17);
        let k18 = Self::from_raw(K18);
        let k19 = Self::from_raw(K19);
        let k20 = Self::from_raw(K20);
        let k21 = Self::from_raw(K21);
        let k22 = Self::from_raw(K22);
        let k23 = Self::from_raw(K23);
        let k24 = Self::from_raw(K24);
        let k25 = Self::from_raw(K25);
        let k26 = Self::from_raw(K26);

        let r = horner!(
            x,
            x,
            [
                k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18, k19, k20, k21,
                k22, k23, k24, k25, k26
            ]
        );
        (r, k1, k2, k3)
    }

    #[inline]
    fn special_poly(x: Self) -> Self {
        // GENERATE: gamma::special_poly f64 22 2.3
        const K0: u64 = 0x3FB5555555555555; // 8.333333333333333e-2
        const K1: u64 = 0x3F6C71C71C71C71D; // 3.4722222222222225e-3
        const K2: u64 = 0xBF65F7268EDAB562; // -2.681327160493894e-3
        const K3: u64 = 0xBF2E13CE46596435; // -2.294720936102682e-4
        const K4: u64 = 0x3F49B0FF67EA3788; // 7.840392207343208e-4
        const K5: u64 = 0x3F12476133E14208; // 6.972819116524847e-5
        const K6: u64 = 0xBF436777EDF04343; // -5.921683876995782e-4
        const K7: u64 = 0xBF0B16B496C91C91; // -5.166758169088918e-5
        const K8: u64 = 0x3F4B7A31BED21A9D; // 8.385413072041385e-4
        const K9: u64 = 0x3F16813B38F7D790; // 8.584903684280402e-5
        const K10: u64 = 0xBF60F0E69E331F78; // -2.067995477407377e-3
        const K11: u64 = 0x3F5332A82C88CF38; // 1.1717455219128602e-3
        const K12: u64 = 0xBF6608165FCD7932; // -2.6894032475147872e-3
        const K13: u64 = 0x3FA90EE5D7B24D31; // 4.894178636564062e-2
        const K14: u64 = 0xBFCD1D6674C1002D; // -2.274597234809347e-1
        const K15: u64 = 0x3FE30A7D96B4D425; // 5.950305884822497e-1
        const K16: u64 = 0xBFF078B777C664BD; // -1.0294718435964534e0
        const K17: u64 = 0x3FF3492E51A7DFF0; // 1.205366438834968e0
        const K18: u64 = 0xBFEC517B561AF98F; // -8.849465066667096e-1
        const K19: u64 = 0x3FD14A384A88C2EB; // 2.701550224183353e-1
        const K20: u64 = 0x3FC344A179080DDD; // 1.5053194436778386e-1
        const K21: u64 = 0xBFC6DF751C10B452; // -1.7869438047765357e-1
        const K22: u64 = 0x3FAB8C6B5443C33F; // 5.3805688892572416e-2

        let k0 = Self::from_raw(K0);
        let k1 = Self::from_raw(K1);
        let k2 = Self::from_raw(K2);
        let k3 = Self::from_raw(K3);
        let k4 = Self::from_raw(K4);
        let k5 = Self::from_raw(K5);
        let k6 = Self::from_raw(K6);
        let k7 = Self::from_raw(K7);
        let k8 = Self::from_raw(K8);
        let k9 = Self::from_raw(K9);
        let k10 = Self::from_raw(K10);
        let k11 = Self::from_raw(K11);
        let k12 = Self::from_raw(K12);
        let k13 = Self::from_raw(K13);
        let k14 = Self::from_raw(K14);
        let k15 = Self::from_raw(K15);
        let k16 = Self::from_raw(K16);
        let k17 = Self::from_raw(K17);
        let k18 = Self::from_raw(K18);
        let k19 = Self::from_raw(K19);
        let k20 = Self::from_raw(K20);
        let k21 = Self::from_raw(K21);
        let k22 = Self::from_raw(K22);

        k0 + horner!(
            x,
            x,
            [
                k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18,
                k19, k20, k21, k22
            ]
        )
    }
}
