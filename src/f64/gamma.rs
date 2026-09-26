use super::f64x2::F64x2;
use super::log_core::log_core_f64x2;
use crate::generic::{reduce_half_revs, round_fi};
use crate::traits::Float as _;

// GENERATE: consts F64x2 PI LN_2 LN_PI
const PI: F64x2 = F64x2::from_bits(0x400921FB54442D18, 0x3CA1A62633145C07); // 3.141592653589793238462643383280e0
const LN_2: F64x2 = F64x2::from_bits(0x3FE62E42FEFA39EF, 0x3C7ABC9E3B39803F); // 6.931471805599453094172321214582e-1
const LN_PI: F64x2 = F64x2::from_bits(0x3FF250D048E7A1BD, 0x3C67ABF2AD8D5088); // 1.144729885849400174143427351353e0

impl crate::generic::Gamma for f64 {
    fn gamma_finite(x: Self) -> Self {
        if x >= 172.0 {
            // Overflow
            return f64::INFINITY;
        } else if x <= -185.0 {
            // Underflow
            let xi = (-x) as u64;
            return 0.0.set_sign((xi & 1) == 0);
        }

        // Split x = k + f + i, such as:
        //  * k is some constant
        //  * i is an integer
        //  * -0.5 <= f <= 0.5
        // so, Γ(x) = Γ(k + f + i)
        //
        // when i = 0:
        //   Γ(x) = Γ(k + f)
        //
        // when i > 0:
        //  Γ(x) = Γ(k + f + i) = Γ(k + f) * prod(j = 0 to i - 1, k + f + j)
        //
        // when i < 0:
        //  Γ(k + f) = Γ(k + f + i) * prod(j = 0 to |i| - 1, k + f + i + j)
        //           = Γ(k + f + i) * prod(j = 0 to |i| - 1, x + j)
        //  Γ(x) = Γ(k + f + i) = Γ(k + f) / prod(j = 0 to |i| - 1, x + j)

        let x = F64x2::new1(x);
        let k = F64x2::new1(2.875);

        let d = x - k;
        let (i_f, i) = round_fi(d.hi());
        let f = d - i_f;

        // gf = Γ(k + f)
        let gf = {
            // GENERATE: gamma_poly F64x2 31 2.875 -0.5 0.5
            const K0: F64x2 = F64x2::from_bits(0x3FFC9A76BE577123, 0xBC858BAA4910955D); // 1.787710898896940310648430613841e0
            const K1: F64x2 = F64x2::from_bits(0x3FF8F2754DDCF90C, 0xBC5E952C6707A0FC); // 1.559193901207950532258122646569e0
            const K2: F64x2 = F64x2::from_bits(0x3FF0D1191949418A, 0xBC97649F30E26B05); // 1.051049326681182812689103810152e0
            const K3: F64x2 = F64x2::from_bits(0x3FDE1F42CF0AE65D, 0x3C58A32EB4779FCE); // 4.706580182933971010639221637751e-1
            const K4: F64x2 = F64x2::from_bits(0x3FC82B358A3AE617, 0xBC6E1A402677FC0B); // 1.888186383201150988955183496664e-1
            const K5: F64x2 = F64x2::from_bits(0x3FAE1F2B30CB1213, 0xBC311B31F344508B); // 5.883154841061268615602355174982e-2
            const K6: F64x2 = F64x2::from_bits(0x3F9240F6D3D6389B, 0x3C228FDEE122F008); // 1.782594364117838206881671768259e-2
            const K7: F64x2 = F64x2::from_bits(0x3F71522CA0D75D02, 0x3C00FD0971ADF4CC); // 4.228758172266868416088314128702e-3
            const K8: F64x2 = F64x2::from_bits(0x3F51FD00825F0768, 0xBBC9FF362C4AAD7A); // 1.097918031050382550820556534275e-3
            const K9: F64x2 = F64x2::from_bits(0x3F298088605FC5DA, 0x3BC7E4281162ACB9); // 1.945654368565159290975867946163e-4
            const K10: F64x2 = F64x2::from_bits(0x3F0B3F4464718C12, 0xBB8E2CA2A67D9F18); // 5.196979014312359439905333492771e-5
            const K11: F64x2 = F64x2::from_bits(0x3ED49E295177A2A9, 0xBB44BB400A1197CC); // 4.915670863671916558241066655308e-6
            const K12: F64x2 = F64x2::from_bits(0x3EC4815453DA85C2, 0x3B613D5FA041FB9C); // 2.444409488003909546465012630043e-6
            const K13: F64x2 = F64x2::from_bits(0xBE84192C940EE4AE, 0xBB0F6ACF4021334E); // -1.497442756717934397544905117778e-7
            const K14: F64x2 = F64x2::from_bits(0x3E86251E6A94B88B, 0x3B254F96FCA7EA27); // 1.649930727958008534274919944674e-7
            const K15: F64x2 = F64x2::from_bits(0xBE65B905492B1DA7, 0x3B07E5B6C8B1ED10); // -4.046175052434777115210896690605e-8
            const K16: F64x2 = F64x2::from_bits(0x3E51CB1491256E6C, 0xBADB1018DDFEF314); // 1.657128574063061176010708832409e-8
            const K17: F64x2 = F64x2::from_bits(0xBE3755CF573D9163, 0x3AA84580E2840DE4); // -5.433148476126354770703322615975e-9
            const K18: F64x2 = F64x2::from_bits(0x3E20A33BF9BA21B3, 0xBABBA9EEF7B67501); // 1.936875505308665426560633680090e-9
            const K19: F64x2 = F64x2::from_bits(0xBE06F9382E90DEB1, 0x3AA0E701BDCE6F19); // -6.686172424155318897092563333409e-10
            const K20: F64x2 = F64x2::from_bits(0x3DF00A1AE2E51F6D, 0xBA9A7239262614F9); // 2.334050478022566203217894712360e-10
            const K21: F64x2 = F64x2::from_bits(0xBDD64DEF09EB7CAB, 0x3A633CF9BADF6B6F); // -8.114303888526165268655674041413e-11
            const K22: F64x2 = F64x2::from_bits(0x3DBF0DB98C25EBCC, 0xBA2A8B77051B6168); // 2.824309602097493615856032550198e-11
            const K23: F64x2 = F64x2::from_bits(0xBDA59B310C5B9E4F, 0x3A4785DFB0BA431C); // -9.825370019485202600965137632928e-12
            const K24: F64x2 = F64x2::from_bits(0x3D8E109223973F55, 0xBA21C00CACADE901); // 3.417964069806284725467197789410e-12
            const K25: F64x2 = F64x2::from_bits(0xBD74EA491D07ED96, 0xBA11C4A35F0CAB9B); // -1.188890230536239848197845584848e-12
            const K26: F64x2 = F64x2::from_bits(0x3D5D2B6C8D7B4027, 0xB9F237BC93C51AFF); // 4.145253048501754407334992833745e-13
            const K27: F64x2 = F64x2::from_bits(0xBD445445FE770C53, 0xB9D0DAC9817EC0DC); // -1.444476042688559496944549897069e-13
            const K28: F64x2 = F64x2::from_bits(0x3D2B4ABC9E3C81C5, 0xB9B734F39E04EA87); // 4.848022531683235880698229322303e-14
            const K29: F64x2 = F64x2::from_bits(0xBD12A8E77F8C8007, 0xB9AEB6A0C2F732D2); // -1.657321602868227395205881671483e-14
            const K30: F64x2 = F64x2::from_bits(0x3D0165E5AA2E9046, 0xB9677BE1A6076E0E); // 7.726279907561184221452103742440e-15
            const K31: F64x2 = F64x2::from_bits(0xBCE98F2CBE5453D6, 0xB98A5FC9AD04353B); // -2.837649724222766200585212298920e-15

            K0 + horner!(
                f,
                f,
                [
                    K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17,
                    K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31
                ]
            )
        };

        let mut descale = 0;

        let y = match i.cmp(&0) {
            core::cmp::Ordering::Equal => gf,
            core::cmp::Ordering::Greater => {
                // gi = prod(j = 0 to i - 1, k + f + j)
                let mut v = k + f;
                let mut gi = v;
                for _ in 1..i {
                    v += 1.0;
                    gi *= v;
                    if gi.hi().abs() > 1e150 {
                        gi = gi.scalbn_medium(-500);
                        descale += 500;
                    }
                }
                // Γ(x) = Γ(k + f) * prod(j = 0 to i - 1, k + f + j) = gf * gi
                gf * gi
            }
            core::cmp::Ordering::Less => {
                // gi = prod(j = 0 to |i| - 1, x + j)
                let mut v = x;
                let mut gi = v;
                if gi.hi().abs() < 1e-150 {
                    gi = gi.scalbn_medium(500);
                    descale += 500;
                }
                for _ in 1..-i {
                    v += 1.0;
                    gi *= v;
                    if gi.hi().abs() > 1e150 {
                        gi = gi.scalbn_medium(-500);
                        descale -= 500;
                    }
                }
                // Γ(x) = Γ(k + f) / prod(j = 0 to |i| - 1, x + j) = gf / gi
                gf / gi
            }
        };

        y.scalbn_to_f64(descale)
    }

    fn ln_gamma_finite(x: Self) -> (Self, i8) {
        fn ln(x: F64x2) -> F64x2 {
            let xm1 = x - 1.0;
            if xm1.hi().exponent() <= -10 {
                // GENERATE: ln_1p_poly F64x2 9 -0.00196 0.00196
                const K2: F64x2 = F64x2::from_bits(0xBFE0000000000000, 0xB9B3ABB47A40BEE4); // -5.000000000000000000000000000010e-1
                const K3: F64x2 = F64x2::from_bits(0x3FD5555555555555, 0x3C7555553A9AEFC5); // 3.333333333333333333333319515223e-1
                const K4: F64x2 = F64x2::from_bits(0xBFD0000000000000, 0x3B18FD3798B2A570); // -2.499999999999999999999948323691e-1
                const K5: F64x2 = F64x2::from_bits(0x3FC999999999999A, 0xBC618AA84C4AAF3C); // 2.000000000000000034947606770695e-1
                const K6: F64x2 = F64x2::from_bits(0xBFC5555555555556, 0x3C6A572B35303A58); // -1.666666666666666737470114920627e-1
                const K7: F64x2 = F64x2::from_bits(0x3FC249249247A1EA, 0x3C459157FFBA48FF); // 1.428571428543954936543353017974e-1
                const K8: F64x2 = F64x2::from_bits(0xBFBFFFFFFFFBB42C, 0x3C5D6CB0144883D8); // -1.249999999960926191955482367882e-1
                const K9: F64x2 = F64x2::from_bits(0x3FBC71D5644B1151, 0x3C57B62EE0DB90AD); // 1.111119623047376515980546472391e-1
                const K10: F64x2 = F64x2::from_bits(0xBFB999A970DADB7C, 0x3C50597C0DD845D7); // -1.000009441876112368316239272785e-1

                let xm1_2 = xm1.square();
                xm1 + horner!(xm1_2, xm1, [K2, K3, K4, K5, K6, K7, K8, K9, K10])
            } else {
                // GENERATE: ln_1p_poly F64x2 11 -0.0079 0.0079
                const K2: F64x2 = F64x2::from_bits(0xBFE0000000000000, 0x39F601E7049BC657); // -4.999999999999999999999999999826e-1
                const K3: F64x2 = F64x2::from_bits(0x3FD5555555555555, 0x3C75555578182061); // 3.333333333333333333333351304303e-1
                const K4: F64x2 = F64x2::from_bits(0xBFD0000000000000, 0xBB23BCD821F56DE1); // -2.500000000000000000000081633102e-1
                const K5: F64x2 = F64x2::from_bits(0x3FC999999999999A, 0xBC6A8A84FF25C652); // 1.999999999999999995918664608235e-1
                const K6: F64x2 = F64x2::from_bits(0xBFC5555555555555, 0xBC6307D829562E3C); // -1.666666666666666656680335369008e-1
                const K7: F64x2 = F64x2::from_bits(0x3FC24924924928CF, 0xBC5E70CFF949F59F); // 1.428571428571729574114667860838e-1
                const K8: F64x2 = F64x2::from_bits(0xBFC000000000073C, 0xBC5C9BF647CF0E65); // -1.250000000000514095296772913210e-1
                const K9: F64x2 = F64x2::from_bits(0x3FBC71C71848932F, 0x3C17D12AFD53B65E); // 1.111111101423147704444889502914e-1
                const K10: F64x2 = F64x2::from_bits(0xBFB9999994168007, 0xBC5A36D09652E250); // -9.999999871661212915701685141442e-2
                const K11: F64x2 = F64x2::from_bits(0x3FB746BEEBD5E64E, 0xBC545A459CB72DBD); // 9.092324502662415655604764984116e-2
                const K12: F64x2 = F64x2::from_bits(0xBFB55657E5167FFB, 0xBC58B6BE62A54AAA); // -8.334874480215142248043222237371e-2

                let (k, lo, ln_hi) = log_core_f64x2(x, 0);
                let lo2 = lo.square();
                let ln_lo = lo + horner!(lo2, lo, [K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12]);

                // ln(x) = ln(2^k * m) = k * ln(2) + ln(m)
                LN_2 * k + (ln_hi + ln_lo)
            }
        }

        fn sin(x: F64x2) -> F64x2 {
            // GENERATE: sin_poly F64x2 9
            const K3: F64x2 = F64x2::from_bits(0xBFC5555555555555, 0xBC65555552CA81AF); // -1.666666666666666666666666009500e-1
            const K5: F64x2 = F64x2::from_bits(0x3F81111111111111, 0x3C0110E7CB9CA23D); // 8.333333333333333333329065980117e-3
            const K7: F64x2 = F64x2::from_bits(0xBF2A01A01A01A01A, 0xBB572C66C6063195); // -1.984126984126984126029916072419e-4
            const K9: F64x2 = F64x2::from_bits(0x3EC71DE3A556C731, 0x3B475F97DC1DE3C8); // 2.755731922398588019213703214996e-6
            const K11: F64x2 = F64x2::from_bits(0xBE5AE64567F53D3B, 0xBAE00414E69D7B82); // -2.505210838543523223314376036023e-8
            const K13: F64x2 = F64x2::from_bits(0x3DE61246139A1624, 0x3A819409726D3F63); // 1.605904383439241982126200763062e-10
            const K15: F64x2 = F64x2::from_bits(0xBD6AE7F3C61EDD83, 0xB9FE33789227A734); // -7.647163171398073276859390202776e-13
            const K17: F64x2 = F64x2::from_bits(0x3CE9529972C79725, 0xB98E162BB6739F98); // 2.811379344450285307153501050905e-15
            const K19: F64x2 = F64x2::from_bits(0xBC62D1434218BDB5, 0xB8F94EABF8C41E64); // -8.160760491214478515701493028147e-18

            let x2 = x * x;
            let x3 = x2 * x;

            x + horner!(x3, x2, [K3, K5, K7, K9, K11, K13, K15, K17, K19])
        }

        fn cos(x: F64x2) -> F64x2 {
            // GENERATE: cos_poly F64x2 8
            const K4: F64x2 = F64x2::from_bits(0x3FA5555555555555, 0x3C4555424437931C); // 4.166666666666666666663512319225e-2
            const K6: F64x2 = F64x2::from_bits(0xBF56C16C16C16C17, 0x3BEFF162A57E0BE3); // -1.388888888888888887780931530771e-3
            const K8: F64x2 = F64x2::from_bits(0x3EFA01A01A01A016, 0xBB9B2646754AC2E7); // 2.480158730158728657597925339703e-5
            const K10: F64x2 = F64x2::from_bits(0xBE927E4FB7789792, 0xBB29D6F71F4B35F6); // -2.755731922397533319638598669759e-7
            const K12: F64x2 = F64x2::from_bits(0x3E21EED8EFE8FA9C, 0xBA9AA526C2C48394); // 2.087675698356730672524610350407e-9
            const K14: F64x2 = F64x2::from_bits(0xBDA9397481DD9560, 0x3A432D3A3990DD4A); // -1.147074454371334149048113104850e-11
            const K16: F64x2 = F64x2::from_bits(0x3D2AE7BB51F3606E, 0xB9CD962E42784608); // 4.779323963823227459005447965865e-14
            const K18: F64x2 = F64x2::from_bits(0xBCA6556752093B7B, 0xB94DA7E65CAE91BA); // -1.549705349810074051809391126314e-16

            let x2 = x * x;
            let x4 = x2 * x2;

            let t = horner!(x4, x2, [K4, K6, K8, K10, K12, K14, K16, K18]);
            (t - x2 * 0.5) + 1.0
        }

        let (y, sign) = if x.abs() <= 45.0 {
            // Split x = k + f + i, such as:
            //  * k is some constant
            //  * i is an integer
            //  * -0.5 <= f <= 0.5
            // so, Γ(x) = Γ(k + f + i)
            //
            // when i = 0:
            //   Γ(x) = Γ(k + f)
            //
            // when i > 0:
            //  Γ(x) = Γ(k + f + i) = Γ(k + f) * prod(j = 0 to i - 1, k + f + j)
            //
            // when i < 0:
            //  Γ(k + f) = Γ(k + f + i) * prod(j = 0 to |i| - 1, k + f + i + j)
            //           = Γ(k + f + i) * prod(j = 0 to |i| - 1, x + j)
            //  Γ(x) = Γ(k + f + i) = Γ(k + f) / prod(j = 0 to |i| - 1, x + j)

            let x = F64x2::new1(x);
            let k = F64x2::new1(2.0);

            let d = x - k;
            let (i_f, i) = round_fi(d.hi());
            let f = d - i_f;

            // lgf = ln(Γ(k + f))
            let lgf = {
                // GENERATE: ln_gamma_poly F64x2 36 2 -0.5 0.50001
                const K1: F64x2 = F64x2::from_bits(0x3FDB0EE6072093CE, 0x3C56CB90701FBFAB); // 4.227843350984671393934879099176e-1
                const K2: F64x2 = F64x2::from_bits(0x3FD4A34CC4A60FA6, 0x3C71873D89122009); // 3.224670334241132182362075833230e-1
                const K3: F64x2 = F64x2::from_bits(0xBFB13E001A557607, 0x3C5FB68BE2F88056); // -6.735230105319809513324605383750e-2
                const K4: F64x2 = F64x2::from_bits(0x3F951322AC7D8483, 0x3C3AFC89088D3161); // 2.058080842778454787900092414132e-2
                const K5: F64x2 = F64x2::from_bits(0xBF7E404FC218F5F2, 0x3C1E4A627D22301F); // -7.385551028673985266273097139096e-3
                const K6: F64x2 = F64x2::from_bits(0x3F67ADD6EADB6C30, 0xBBF5B78293ECFA5F); // 2.890510330741523285752986798756e-3
                const K7: F64x2 = F64x2::from_bits(0xBF538AC5C2BF8E08, 0x3BE8A4C0CD1FD7D2); // -1.192753911703260977113961817521e-3
                const K8: F64x2 = F64x2::from_bits(0x3F40B36AF86396E9, 0xBBE069869F520F8E); // 5.096695247430424223358306313901e-4
                const K9: F64x2 = F64x2::from_bits(0xBF2D3FD4C76D2FC8, 0x3BBC7F30637F98AD); // -2.231547584535793797590581821663e-4
                const K10: F64x2 = F64x2::from_bits(0x3F1A127B0F17D65A, 0x3BA9B6784172BC37); // 9.945751278180853370278043874578e-5
                const K11: F64x2 = F64x2::from_bits(0xBF078DE5BD7C81EF, 0x3B70548F1EEBD550); // -4.492623673813314182985614684945e-5
                const K12: F64x2 = F64x2::from_bits(0x3EF580DCEE66EB02, 0x3B92C4F616A00980); // 2.050721277567069206039175050237e-5
                const K13: F64x2 = F64x2::from_bits(0xBEE3CBC963CE2240, 0x3B80874D36C62932); // -9.439488275268391195481996641389e-6
                const K14: F64x2 = F64x2::from_bits(0x3ED2597A39F34A9A, 0x3B05EDEDC0CB9EC7); // 4.374866789907472930074701780064e-6
                const K15: F64x2 = F64x2::from_bits(0xBEC11B2EB767965B, 0x3B49E39F79E2CE51); // -2.039215753801485578508317068122e-6
                const K16: F64x2 = F64x2::from_bits(0x3EB0064CDEB234CE, 0xBB59DCBE3946131D); // 9.551412130410533460958461155424e-7
                const K17: F64x2 = F64x2::from_bits(0xBE9E2600D93C5BBE, 0x3B3AAD572383E8C2); // -4.492469198742686574396636298584e-7
                const K18: F64x2 = F64x2::from_bits(0x3E8C76BBB3EDB71A, 0xBB284AC2C6B3DFC1); // 2.120718480507544680937399579938e-7
                const K19: F64x2 = F64x2::from_bits(0xBE7AF5A6CBE1CB09, 0xBB169E3877737C16); // -1.004322482693896438213232902979e-7
                const K20: F64x2 = F64x2::from_bits(0x3E699B93C28688CC, 0xBB08C8406C21D5B8); // 4.769810174893028071660751379732e-8
                const K21: F64x2 = F64x2::from_bits(0xBE5862C72F692EEF, 0x3ACCE93DE6C5BB05); // -2.271109430578690504772161414684e-8
                const K22: F64x2 = F64x2::from_bits(0x3E47469D9B61F70F, 0x3AE569F4185DAB32); // 1.083865873115697167381853909414e-8
                const K23: F64x2 = F64x2::from_bits(0xBE36434B2CFAA05D, 0x3AD442BB0AE4723B); // -5.183477383138976749176573867046e-9
                const K24: F64x2 = F64x2::from_bits(0x3E2555AA48ABE539, 0xBAC4129084604839); // 2.483677768113273727948991063949e-9
                const K25: F64x2 = F64x2::from_bits(0xBE147B07123F8A42, 0x3A780F4B56FAC2C0); // -1.192126460908977237512502955673e-9
                const K26: F64x2 = F64x2::from_bits(0x3E03B138614E04CA, 0x3AA0DDC25F6E068E); // 5.731203867472705309374566811387e-10
                const K27: F64x2 = F64x2::from_bits(0xBDF2F7A8BC265BEF, 0x3A680D285B5B4A9A); // -2.760122652568227958083301686425e-10
                const K28: F64x2 = F64x2::from_bits(0x3DE24B63D0396DD7, 0xBA81BBD2A8D93119); // 1.331099467734351145993975832208e-10
                const K29: F64x2 = F64x2::from_bits(0xBDD19A226B1B342F, 0x3A7CE071A589647D); // -6.403602195081112747904879736530e-11
                const K30: F64x2 = F64x2::from_bits(0x3DC0F86EE2AA5CA7, 0x3A0EBBA7DC467B63); // 3.086905412985665043714690981573e-11
                const K31: F64x2 = F64x2::from_bits(0xBDB10138281FCB26, 0xBA50DF44EB88D91E); // -1.546574197507270976808910881965e-11
                const K32: F64x2 = F64x2::from_bits(0x3DA0C515DBA7298B, 0xBA4C48A171B5ECE7); // 7.626051582023145483915710453512e-12
                const K33: F64x2 = F64x2::from_bits(0xBD890EE290F302EF, 0xBA118BEAC09B56F9); // -2.848781221698287440532080399828e-12
                const K34: F64x2 = F64x2::from_bits(0x3D75F2B4F03A237C, 0xBA09EA007B2D5C03); // 1.247603529502751232341366700243e-12
                const K35: F64x2 = F64x2::from_bits(0xBD79330AEAD7411A, 0x3A16B817A97EF17A); // -1.432419215660462389905381982977e-12
                const K36: F64x2 = F64x2::from_bits(0x3D6AAB78B8436D26, 0x39D09B400E90651D); // 7.580016127701884965269890955627e-13

                horner!(
                    f,
                    f,
                    [
                        K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17,
                        K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32,
                        K33, K34, K35, K36
                    ]
                )
            };

            match i.cmp(&0) {
                core::cmp::Ordering::Equal => (lgf, 1),
                core::cmp::Ordering::Greater => {
                    // gi = prod(j = 0 to i - 1, k + f + j)
                    let mut descale = 0i16;
                    let mut v = k + f;
                    let mut gi = v;
                    for _ in 1..i {
                        v += 1.0;
                        gi *= v;
                        if gi.hi().abs() > 1e150 {
                            gi = gi.scalbn_medium(-500);
                            descale += 500;
                        }
                    }
                    // ln(abs(Γ(x))) = ln(Γ(k + f)) + ln(prod(j = 0 to i - 1, k + f + j)) = lgf + ln(gi)
                    (lgf + ln(gi) + LN_2 * f64::from(descale), 1)
                }
                core::cmp::Ordering::Less => {
                    // gi = prod(j = 0 to |i| - 1, x + j)
                    let mut descale = 0i16;
                    let mut v = x;
                    let mut gi = v;
                    if gi.hi().abs() < 1e-150 {
                        gi = gi.scalbn_medium(500);
                        descale += 500;
                    }
                    for _ in 1..-i {
                        v += 1.0;
                        gi *= v;
                        if gi.hi().abs() > 1e150 {
                            gi = gi.scalbn_medium(-500);
                            descale -= 500;
                        }
                    }
                    // ln(abs(Γ(x))) = ln(Γ(k + f)) - ln(abs(prod(j = 0 to |i| - 1, x + j))) = lgf - ln(abs(gi))
                    let sign = if gi.hi().is_sign_negative() { -1 } else { 1 };
                    (lgf - ln(gi.abs()) + LN_2 * f64::from(descale), sign)
                }
            }
        } else {
            // Use Lanczos approximation:
            // Γ(x) = √(2π) * (x + g - 0.5)^(x - 0.5) * exp(-(x + g - 0.5)) * Ag(x - 1)
            //      = (x + g - 0.5)^(x - 0.5) * exp(-(x + g - 0.5)) * P(1 / x)
            // with let P(1 / x) = √(2π) * Ag(x - 1) = Γ(x) / ((x + g - 0.5)^(x - 0.5) * exp(-(x + g - 0.5)))
            //
            // choose g = 0.5, so:
            // ln(Γ(x)) = (x - 0.5) * ln(x) - x + ln(P(1 / x))
            //
            // For x < 0.5, use reflection formula:
            // Γ(x)*Γ(1-x) = π/sin(πx) => Γ(x) = π/(sin(πx)*Γ(1-x))

            // nx = x or 1 - x, so nx >= 0.5
            let reflect = x.is_sign_negative();
            let nx = if reflect {
                F64x2::sub11(1.0, x)
            } else {
                F64x2::new1(x)
            };

            // p = P(1 / nx)
            let p = {
                // GENERATE: gamma_lanczos_poly F64x2 9 0.5 2.93e-39 0.022223
                const K0: F64x2 = F64x2::from_bits(0x40040D931FF62706, 0xBCAA6A0D6F734365); // 2.506628274631000502415765307447e0
                const K1: F64x2 = F64x2::from_bits(0x3FCABCC42A9D895D, 0xBC5BC5B509658A44); // 2.088856895525833752011156633309e-1
                const K2: F64x2 = F64x2::from_bits(0x3F81D32D71BE5B93, 0x3C269ED00C83FB3C); // 8.703570398024307585854645818948e-3
                const K3: F64x2 = F64x2::from_bits(0xBF7B8792FC787D15, 0xBC1B63009BB23951); // -6.721090474030041162568244073157e-3
                const K4: F64x2 = F64x2::from_bits(0xBF42D92340C82F32, 0x3BD75859BD27AD1E); // -5.752012380652925138697500384312e-4
                const K5: F64x2 = F64x2::from_bits(0x3F601985A6CD5224, 0x3BF4A8C31E99B8DE); // 1.965294874404459741276072510530e-3
                const K6: F64x2 = F64x2::from_bits(0x3F26E8C0502F1D3B, 0xBBC2E21A771EEB40); // 1.747832066867157454773808239276e-4
                const K7: F64x2 = F64x2::from_bits(0xBF5851F0FA6CDAE2, 0x3BEFAC1912E012A6); // -1.484380083366758742765523474202e-3
                const K8: F64x2 = F64x2::from_bits(0xBF20D52DE91D4F7B, 0xBBA57C94D79C2952); // -1.284235518662599273329735872778e-4
                const K9: F64x2 = F64x2::from_bits(0x3F611F449B9249DD, 0x3BC3637631223E74); // 2.090105056205644835111295019112e-3

                let inv_nx = nx.recip();
                K0 + horner!(inv_nx, inv_nx, [K1, K2, K3, K4, K5, K6, K7, K8, K9])
            };

            // ln(Γ(nx)) = (nx - 0.5) * ln(nx) - nx + ln(P(1 / nx))
            //           = nx * (ln(nx) - 1) - 0.5 * ln(nx) + ln(P(1 / nx))
            // The second form is used to avoid overflow in intermediate operations.
            let ln_nx = ln(nx);
            let lg_nx = (ln_nx - 1.0) * nx - ln_nx.halve() + ln(p);

            if reflect {
                let (n, z) = reduce_half_revs(x);
                let z_rad = PI * z;
                let sinpix = match n {
                    0 => sin(z_rad),
                    1 => cos(z_rad),
                    2 => -sin(z_rad),
                    3 => -cos(z_rad),
                    _ => unreachable!(),
                };

                // ln(abs(Γ(x))) = ln(π) - ln(abs(sin(πx))) - ln(Γ(1-x))
                let lgx = LN_PI - ln(sinpix.abs()) - lg_nx;
                let sign = if sinpix.hi().is_sign_negative() {
                    -1
                } else {
                    1
                };

                (lgx, sign)
            } else {
                // ln(abs(Γ(x))) = ln(Γ(nx))
                (lg_nx, 1)
            }
        };

        // F64x2 can produce NaN on overflow
        if !y.hi().is_finite() {
            (f64::INFINITY, sign)
        } else {
            (y.to_f64(), sign)
        }
    }
}
