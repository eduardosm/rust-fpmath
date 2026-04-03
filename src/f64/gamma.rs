use crate::double::NormDouble;

impl crate::generic::Gamma for f64 {
    #[inline]
    fn lo_th() -> Self {
        -10000.0
    }

    #[inline]
    fn hi_th() -> Self {
        10000.0
    }

    #[inline]
    fn th_1() -> Self {
        1.1
    }

    #[inline]
    fn th_2() -> Self {
        2.3
    }

    #[inline]
    fn th_3() -> Self {
        7.0
    }

    const POLY_OFF: u8 = 5;

    #[inline]
    fn half_ln_2_pi() -> NormDouble<Self> {
        // GENERATE: gamma::consts f64
        const HALF_LN_2_PI_HI: f64 = f64::from_bits(0x3FED67F1C864BEB4); // 9.189385332046727e-1
        const HALF_LN_2_PI_LO: f64 = f64::from_bits(0x3C94D252F2400510); // 7.223936088184323e-17

        NormDouble::with_parts(HALF_LN_2_PI_HI, HALF_LN_2_PI_LO)
    }

    #[inline]
    fn lgamma_poly_1(x: Self) -> (Self, Self, Self, Self) {
        // GENERATE: gamma::lgamma_poly f64 26 1 0.5 1.1001
        const K1: f64 = f64::from_bits(0xBFE2788CFC6FB619); // -5.772156649015329e-1
        const K2: f64 = f64::from_bits(0x3FEA51A6625307D3); // 8.224670334241132e-1
        const K3: f64 = f64::from_bits(0xBFD9A4D55BEAB2D5); // -4.006856343865313e-1
        const K4: f64 = f64::from_bits(0x3FD151322AC7D84A); // 2.7058080842778465e-1
        const K5: f64 = f64::from_bits(0xBFCA8B9C17AA7467); // -2.0738555102880982e-1
        const K6: f64 = f64::from_bits(0x3FC5B40CB1005A28); // 1.6955717699666306e-1
        const K7: f64 = f64::from_bits(0xBFC2703A1DB0425A); // -1.4404989671358487e-1
        const K8: f64 = f64::from_bits(0x3FC010B36C29B9A6); // 1.255096700801463e-1
        const K9: f64 = f64::from_bits(0xBFBC80672DDE4D26); // -1.1133427495797274e-1
        const K10: f64 = f64::from_bits(0x3FB9A01BA4EC8538); // 1.0009930397474853e-1
        const K11: f64 = f64::from_bits(0xBFB748BF8D505593); // -9.095380020183645e-2
        const K12: f64 = f64::from_bits(0x3FB557D80D3B0BAA); // 8.337164233220898e-2
        const K13: f64 = f64::from_bits(0xBFB3AB2A85D1DC7A); // -7.683053748211136e-2
        const K14: f64 = f64::from_bits(0x3FB225FF02F2284D); // 7.08922750656324e-2
        const K15: f64 = f64::from_bits(0xBFB384C239C83F7F); // -7.624448691080764e-2
        const K16: f64 = f64::from_bits(0x3F95419F7F53E320); // 2.0758144518342125e-2
        const K17: f64 = f64::from_bits(0x3F9168340914CC99); // 1.699906640460522e-2
        const K18: f64 = f64::from_bits(0x3FFEDBF90E5FE0AA); // 1.928704315329052e0
        const K19: f64 = f64::from_bits(0x4026C176611EDB2C); // 1.1377856287972087e1
        const K20: f64 = f64::from_bits(0x4045190442C39B12); // 4.219544252921163e1
        const K21: f64 = f64::from_bits(0x405A9D2F7FF1BA92); // 1.0645602415663078e2
        const K22: f64 = f64::from_bits(0x4067CD207F00C1BF); // 1.9041021680972878e2
        const K23: f64 = f64::from_bits(0x406DDABAAD9CD985); // 2.3883528786310868e2
        const K24: f64 = f64::from_bits(0x40693ADFB5F6BF8E); // 2.0183980844681759e2
        const K25: f64 = f64::from_bits(0x4059F58E0C89743E); // 1.0383679498122453e2
        const K26: f64 = f64::from_bits(0x4038F16511CF45FF); // 2.4942948449233878e1

        let r = horner!(
            x,
            x,
            [
                K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21,
                K22, K23, K24, K25, K26
            ]
        );
        (r, K1, K2, K3)
    }

    #[inline]
    fn lgamma_poly_2(x: Self) -> (Self, Self, Self, Self) {
        // GENERATE: gamma::lgamma_poly f64 26 2 1.0999 2.3001
        const K1: f64 = f64::from_bits(0x3FDB0EE6072093CE); // 4.2278433509846713e-1
        const K2: f64 = f64::from_bits(0x3FD4A34CC4A60FA6); // 3.224670334241132e-1
        const K3: f64 = f64::from_bits(0xBFB13E001A557607); // -6.73523010531981e-2
        const K4: f64 = f64::from_bits(0x3F951322AC7D8469); // 2.0580808427784456e-2
        const K5: f64 = f64::from_bits(0xBF7E404FC218F52F); // -7.3855510286738165e-3
        const K6: f64 = f64::from_bits(0x3F67ADD6EADBCE6F); // 2.890510330752431e-3
        const K7: f64 = f64::from_bits(0xBF538AC5C2BF075B); // -1.192753911695785e-3
        const K8: f64 = f64::from_bits(0x3F40B36AF811CBF9); // 5.096695241618702e-4
        const K9: f64 = f64::from_bits(0xBF2D3FD4CA9756B1); // -2.2315475989269842e-4
        const K10: f64 = f64::from_bits(0x3F1A127B510D64FD); // 9.945752777919602e-5
        const K11: f64 = f64::from_bits(0xBF078DE385E76226); // -4.492617221147172e-5
        const K12: f64 = f64::from_bits(0x3EF580D235E1E8CA); // 2.050705676779548e-5
        const K13: f64 = f64::from_bits(0xBEE3CC7D4FD7893E); // -9.440797380248623e-6
        const K14: f64 = f64::from_bits(0x3ED258B5A6842E9C); // 4.374151650897832e-6
        const K15: f64 = f64::from_bits(0xBEC1016D13BA0E79); // -2.0272219936628015e-6
        const K16: f64 = f64::from_bits(0x3EB08A3BDBAA3E36); // 9.858593858247807e-7
        const K17: f64 = f64::from_bits(0xBE9F88691A3B69A9); // -4.698761364281987e-7
        const K18: f64 = f64::from_bits(0xBE41B224E84CF21E); // -8.240284241083623e-9
        const K19: f64 = f64::from_bits(0xBE9EC6BD7BBC2C6A); // -4.586030361997586e-7
        const K20: f64 = f64::from_bits(0x3E86034DDEBEEAB6); // 1.640089374291915e-7
        const K21: f64 = f64::from_bits(0x3EB6B1EC9F3AD0E3); // 1.3527284154326575e-6
        const K22: f64 = f64::from_bits(0x3EC59F32DBC06B30); // 2.577527736147797e-6
        const K23: f64 = f64::from_bits(0x3EC5F104E397AA2A); // 2.615628343906737e-6
        const K24: f64 = f64::from_bits(0x3EBB389DFEF83D54); // 1.6225076214826323e-6
        const K25: f64 = f64::from_bits(0x3EA35407E7E24522); // 5.760266075396234e-7
        const K26: f64 = f64::from_bits(0x3E7938471B002083); // 9.395120659239e-8

        let r = horner!(
            x,
            x,
            [
                K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21,
                K22, K23, K24, K25, K26
            ]
        );
        (r, K1, K2, K3)
    }

    #[inline]
    fn special_poly(x: Self) -> Self {
        // GENERATE: gamma::special_poly f64 22 2.3
        const K0: f64 = f64::from_bits(0x3FB5555555555555); // 8.333333333333333e-2
        const K1: f64 = f64::from_bits(0x3F6C71C71C71C71D); // 3.4722222222222225e-3
        const K2: f64 = f64::from_bits(0xBF65F7268EDAB561); // -2.6813271604938936e-3
        const K3: f64 = f64::from_bits(0xBF2E13CE46596AE7); // -2.2947209361031466e-4
        const K4: f64 = f64::from_bits(0x3F49B0FF67EAC1E3); // 7.84039220738161e-4
        const K5: f64 = f64::from_bits(0x3F1247613306CF96); // 6.972819097122852e-5
        const K6: f64 = f64::from_bits(0xBF436777EA57F79C); // -5.921683811604852e-4
        const K7: f64 = f64::from_bits(0xBF0B16B9EFF70F16); // -5.166773734857474e-5
        const K8: f64 = f64::from_bits(0x3F4B7A37979321B7); // 8.385440297838308e-4
        const K9: f64 = f64::from_bits(0x3F167ED0E7CAF987); // 8.58130460512416e-5
        const K10: f64 = f64::from_bits(0xBF60F0219D0C087A); // -2.0676285279260216e-3
        const K11: f64 = f64::from_bits(0x3F5326611B6857FC); // 1.1688183127538379e-3
        const K12: f64 = f64::from_bits(0xBF65E16771C7E50C); // -2.670957603410634e-3
        const K13: f64 = f64::from_bits(0x3FA902CC210FE3C8); // 4.884946731176104e-2
        const K14: f64 = f64::from_bits(0xBFCD115B3FE3439E); // -2.2709217662352627e-1
        const K15: f64 = f64::from_bits(0x3FE300FAA82E5E86); // 5.938695225513164e-1
        const K16: f64 = f64::from_bits(0xBFF06CE1D839B6BB); // -1.026582569723742e0
        const K17: f64 = f64::from_bits(0x3FF3324939542702); // 1.1997768630412229e0
        const K18: f64 = f64::from_bits(0xBFEC0E16BDFF841F); // -8.76719828691538e-1
        const K19: f64 = f64::from_bits(0x3FD0B87FD49A8668); // 2.612609459120905e-1
        const K20: f64 = f64::from_bits(0x3FC41EBB55E86478); // 1.5718785948482483e-1
        const K21: f64 = f64::from_bits(0xBFC74458534C3DE1); // -1.8177322451855019e-1
        const K22: f64 = f64::from_bits(0x3FABE35113B97418); // 5.446866384293275e-2

        K0 + horner!(
            x,
            x,
            [
                K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18,
                K19, K20, K21, K22
            ]
        )
    }
}
