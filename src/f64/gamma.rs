use super::log_core::log_core;
use crate::double::Double;
use crate::generic::round_as_i_f;
use crate::scalbn;
use crate::traits::Float as _;

// GENERATE: consts Double<f64> PI LN_2 LN_PI
const PI: Double<f64> = Double::new(
    f64::from_bits(0x400921FB54442D18),
    f64::from_bits(0x3CA1A62633145C07),
); // 3.141592653589793238462643383280e0
const LN_2: Double<f64> = Double::new(
    f64::from_bits(0x3FE62E42FEFA39EF),
    f64::from_bits(0x3C7ABC9E3B39803F),
); // 6.931471805599453094172321214582e-1
const LN_PI: Double<f64> = Double::new(
    f64::from_bits(0x3FF250D048E7A1BD),
    f64::from_bits(0x3C67ABF2AD8D5088),
); // 1.144729885849400174143427351353e0

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

        let x = Double::from(x);
        let k = Double::from(2.875);

        let d = x - k;
        let (i, i_f) = round_as_i_f(d.hi());
        let f = (d - i_f).normalize();

        // gf = Γ(k + f)
        let gf = {
            // GENERATE: gamma_poly Double<f64> 31 2.875 -0.5 0.5
            const K0: Double<f64> = Double::new(
                f64::from_bits(0x3FFC9A76BE577122),
                f64::from_bits(0x3CAA9D156DBBDAA9),
            ); // 1.787710898896940310648430613841e0
            const K1: Double<f64> = Double::new(
                f64::from_bits(0x3FF8F2754DDCF90B),
                f64::from_bits(0x3CAF0B569CC7C2F8),
            ); // 1.559193901207950532258122646569e0
            const K2: Double<f64> = Double::new(
                f64::from_bits(0x3FF0D11919494189),
                f64::from_bits(0x3CA44DB0678ECA7D),
            ); // 1.051049326681182812689103810152e0
            const K3: Double<f64> = Double::new(
                f64::from_bits(0x3FDE1F42CF0AE65D),
                f64::from_bits(0x3C58A32EB4779FCE),
            ); // 4.706580182933971010639221637751e-1
            const K4: Double<f64> = Double::new(
                f64::from_bits(0x3FC82B358A3AE616),
                f64::from_bits(0x3C70F2DFECC401FA),
            ); // 1.888186383201150988955183496664e-1
            const K5: Double<f64> = Double::new(
                f64::from_bits(0x3FAE1F2B30CB1212),
                f64::from_bits(0x3C5BB933832EEBDD),
            ); // 5.883154841061268615602355174982e-2
            const K6: Double<f64> = Double::new(
                f64::from_bits(0x3F9240F6D3D6389B),
                f64::from_bits(0x3C228FDEE122F008),
            ); // 1.782594364117838206881671768259e-2
            const K7: Double<f64> = Double::new(
                f64::from_bits(0x3F71522CA0D75D02),
                f64::from_bits(0x3C00FD0971ADF4CC),
            ); // 4.228758172266868416088314128702e-3
            const K8: Double<f64> = Double::new(
                f64::from_bits(0x3F51FD00825F0767),
                f64::from_bits(0x3C0E600C9D3B5528),
            ); // 1.097918031050382550820556534275e-3
            const K9: Double<f64> = Double::new(
                f64::from_bits(0x3F298088605FC5DA),
                f64::from_bits(0x3BC7E4281162ACB9),
            ); // 1.945654368565159290975867946163e-4
            const K10: Double<f64> = Double::new(
                f64::from_bits(0x3F0B3F4464718C11),
                f64::from_bits(0x3BBC3A6BAB304C1D),
            ); // 5.196979014312359439905333492771e-5
            const K11: Double<f64> = Double::new(
                f64::from_bits(0x3ED49E295177A2A8),
                f64::from_bits(0x3B8EB44BFF5EE683),
            ); // 4.915670863671916558241066655308e-6
            const K12: Double<f64> = Double::new(
                f64::from_bits(0x3EC4815453DA85C2),
                f64::from_bits(0x3B613D5FA041FB9C),
            ); // 2.444409488003909546465012630043e-6
            const K13: Double<f64> = Double::new(
                f64::from_bits(0xBE84192C940EE4AE),
                f64::from_bits(0xBB0F6ACF4021334E),
            ); // -1.497442756717934397544905117778e-7
            const K14: Double<f64> = Double::new(
                f64::from_bits(0x3E86251E6A94B88B),
                f64::from_bits(0x3B254F96FCA7EA27),
            ); // 1.649930727958008534274919944674e-7
            const K15: Double<f64> = Double::new(
                f64::from_bits(0xBE65B905492B1DA6),
                f64::from_bits(0xBB140D249BA70978),
            ); // -4.046175052434777115210896690605e-8
            const K16: Double<f64> = Double::new(
                f64::from_bits(0x3E51CB1491256E6B),
                f64::from_bits(0x3B0C9DFCE440219E),
            ); // 1.657128574063061176010708832409e-8
            const K17: Double<f64> = Double::new(
                f64::from_bits(0xBE3755CF573D9162),
                f64::from_bits(0xBAEE7BA7F1D7BF22),
            ); // -5.433148476126354770703322615975e-9
            const K18: Double<f64> = Double::new(
                f64::from_bits(0x3E20A33BF9BA21B2),
                f64::from_bits(0x3AD91584421262C0),
            ); // 1.936875505308665426560633680090e-9
            const K19: Double<f64> = Double::new(
                f64::from_bits(0xBE06F9382E90DEB0),
                f64::from_bits(0xBAB78C7F2118C873),
            ); // -6.686172424155318897092563333409e-10
            const K20: Double<f64> = Double::new(
                f64::from_bits(0x3DF00A1AE2E51F6C),
                f64::from_bits(0x3AA2C6E36CECF584),
            ); // 2.334050478022566203217894712360e-10
            const K21: Double<f64> = Double::new(
                f64::from_bits(0xBDD64DEF09EB7CAA),
                f64::from_bits(0xBA8B30C191482524),
            ); // -8.114303888526165268655674041413e-11
            const K22: Double<f64> = Double::new(
                f64::from_bits(0x3DBF0DB98C25EBCB),
                f64::from_bits(0x3A6E57488FAE49EA),
            ); // 2.824309602097493615856032550198e-11
            const K23: Double<f64> = Double::new(
                f64::from_bits(0xBDA59B310C5B9E4E),
                f64::from_bits(0xBA543D1027A2DE72),
            ); // -9.825370019485202600965137632928e-12
            const K24: Double<f64> = Double::new(
                f64::from_bits(0x3D8E109223973F54),
                f64::from_bits(0x3A371FF9A9A90B80),
            ); // 3.417964069806284725467197789410e-12
            const K25: Double<f64> = Double::new(
                f64::from_bits(0xBD74EA491D07ED96),
                f64::from_bits(0xBA11C4A35F0CAB9B),
            ); // -1.188890230536239848197845584848e-12
            const K26: Double<f64> = Double::new(
                f64::from_bits(0x3D5D2B6C8D7B4026),
                f64::from_bits(0x3A06E421B61D7281),
            ); // 4.145253048501754407334992833745e-13
            const K27: Double<f64> = Double::new(
                f64::from_bits(0xBD445445FE770C53),
                f64::from_bits(0xB9D0DAC9817EC0DC),
            ); // -1.444476042688559496944549897069e-13
            const K28: Double<f64> = Double::new(
                f64::from_bits(0x3D2B4ABC9E3C81C4),
                f64::from_bits(0x39DA32C3187EC55E),
            ); // 4.848022531683235880698229322303e-14
            const K29: Double<f64> = Double::new(
                f64::from_bits(0xBD12A8E77F8C8007),
                f64::from_bits(0xB9AEB6A0C2F732D2),
            ); // -1.657321602868227395205881671483e-14
            const K30: Double<f64> = Double::new(
                f64::from_bits(0x3D0165E5AA2E9045),
                f64::from_bits(0x39BF4420F2CFC490),
            ); // 7.726279907561184221452103742440e-15
            const K31: Double<f64> = Double::new(
                f64::from_bits(0xBCE98F2CBE5453D6),
                f64::from_bits(0xB98A5FC9AD04353B),
            ); // -2.837649724222766200585212298920e-15

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
                    v = v + 1.0;
                    gi = (gi * v).normalize();
                    if gi.hi().abs() > 1e150 {
                        gi = gi.pmul1(f64::exp2i_fast(-500));
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
                    gi = gi.pmul1(f64::exp2i_fast(500));
                    descale += 500;
                }
                for _ in 1..-i {
                    v = v + 1.0;
                    gi = (gi * v).normalize();
                    if gi.hi().abs() > 1e150 {
                        gi = gi.pmul1(f64::exp2i_fast(-500));
                        descale -= 500;
                    }
                }
                // Γ(x) = Γ(k + f) / prod(j = 0 to |i| - 1, x + j) = gf / gi
                gf / gi
            }
        };

        scalbn(y.to_single(), descale)
    }

    fn ln_gamma_finite(x: Self) -> (Self, i8) {
        fn ln(x: Double<f64>) -> Double<f64> {
            let xm1 = (x - 1.0).normalize();
            if xm1.hi().abs() < 0.0002 {
                // GENERATE: ln_1p_poly Double<f64> 10 -0.00021 0.00021
                const K2: Double<f64> = Double::new(
                    f64::from_bits(0xBFE0000000000000),
                    f64::from_bits(0xB7AA03095EEEDDFE),
                ); // -5.000000000000000000000000000000e-1
                const K3: Double<f64> = Double::new(
                    f64::from_bits(0x3FD5555555555555),
                    f64::from_bits(0x3C75555555555555),
                ); // 3.333333333333333333333333333333e-1
                const K4: Double<f64> = Double::new(
                    f64::from_bits(0xBFCFFFFFFFFFFFFF),
                    f64::from_bits(0xBC7FFFFFFFFFFFEA),
                ); // -2.499999999999999999999999999999e-1
                const K5: Double<f64> = Double::new(
                    f64::from_bits(0x3FC9999999999999),
                    f64::from_bits(0x3C733333333332FE),
                ); // 1.999999999999999999999999999998e-1
                const K6: Double<f64> = Double::new(
                    f64::from_bits(0xBFC5555555555555),
                    f64::from_bits(0xBC655556A1FB35DD),
                ); // -1.666666666666666666666752654148e-1
                const K7: Double<f64> = Double::new(
                    f64::from_bits(0x3FC2492492492492),
                    f64::from_bits(0x3C624926C9D21018),
                ); // 1.428571428571428571428718132890e-1
                const K8: Double<f64> = Double::new(
                    f64::from_bits(0xBFBFFFFFFFFFFFDF),
                    f64::from_bits(0xBC6C54BAC3217EAD),
                ); // -1.249999999999995543196040978647e-1
                const K9: Double<f64> = Double::new(
                    f64::from_bits(0x3FBC71C71C71C6F1),
                    f64::from_bits(0x3C6EDD6F83F84C7C),
                ); // 1.111111111111105215838822629663e-1
                const K10: Double<f64> = Double::new(
                    f64::from_bits(0xBFB99999C5017CB1),
                    f64::from_bits(0xBC610FF2C8D89D83),
                ); // -1.000000101062024295439204427796e-1
                const K11: Double<f64> = Double::new(
                    f64::from_bits(0x3FB745D1A37352B5),
                    f64::from_bits(0x3C578A135F2C2133),
                ); // 9.090910187235120243080075052703e-2

                let xm1_2 = xm1.square();
                xm1 + horner!(xm1_2, xm1, [K2, K3, K4, K5, K6, K7, K8, K9, K10, K11])
            } else {
                // GENERATE: ln_1p_poly Double<f64> 12 -0.0001 0.015625
                const K2: Double<f64> = Double::new(
                    f64::from_bits(0xBFDFFFFFFFFFFFFF),
                    f64::from_bits(0xBC8FFFFFFFFFE50B),
                ); // -4.999999999999999999999999999575e-1
                const K3: Double<f64> = Double::new(
                    f64::from_bits(0x3FD5555555555555),
                    f64::from_bits(0x3C755555503E169D),
                ); // 3.333333333333333333333330701451e-1
                const K4: Double<f64> = Double::new(
                    f64::from_bits(0xBFCFFFFFFFFFFFFF),
                    f64::from_bits(0xBC7FFFDCF26BB33A),
                ); // -2.499999999999999999995360768238e-1
                const K5: Double<f64> = Double::new(
                    f64::from_bits(0x3FC9999999999999),
                    f64::from_bits(0x3C72C1CF308A5DA0),
                ); // 1.999999999999999996158174813676e-1
                const K6: Double<f64> = Double::new(
                    f64::from_bits(0xBFC555555555554E),
                    f64::from_bits(0xBC79FA6E7B0F7597),
                ); // -1.666666666666664856583174156367e-1
                const K7: Double<f64> = Double::new(
                    f64::from_bits(0x3FC2492492491D0C),
                    f64::from_bits(0x3C7CA4F717D680A3),
                ); // 1.428571428570894168191095072070e-1
                const K8: Double<f64> = Double::new(
                    f64::from_bits(0xBFBFFFFFFFF48A74),
                    f64::from_bits(0xBC3EB7D76C482152),
                ); // -1.249999999895779495550182465243e-1
                const K9: Double<f64> = Double::new(
                    f64::from_bits(0x3FBC71C71685F133),
                    f64::from_bits(0x3C6704E4C74FDD78),
                ); // 1.111111097324666802020255908172e-1
                const K10: Double<f64> = Double::new(
                    f64::from_bits(0xBFB9999782F1106E),
                    f64::from_bits(0xBC60F6BD78D1CBC4),
                ); // -9.999987551515385944711006490894e-2
                const K11: Double<f64> = Double::new(
                    f64::from_bits(0x3FB7455289CDAA2F),
                    f64::from_bits(0x3C42A54791E45EB4),
                ); // 9.090152611103329207782503244107e-2
                const K12: Double<f64> = Double::new(
                    f64::from_bits(0xBFB541E6981AA17C),
                    f64::from_bits(0xBC2B1F268BB2B12B),
                ); // -8.303681577183924596113943880234e-2
                const K13: Double<f64> = Double::new(
                    f64::from_bits(0x3FB1F2FD99832E3A),
                    f64::from_bits(0x3C65614552AECCFF),
                ); // 7.011399267064969781543087974875e-2

                let (k, lo, ln_hi) = log_core(x);
                let lo2 = lo * lo;
                let ln_lo = lo
                    + horner!(
                        lo2,
                        lo,
                        [K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13]
                    );

                // ln(x) = ln(2^k * m) = k * ln(2) + ln(m)
                LN_2 * k + (ln_hi + ln_lo)
            }
        }

        fn sin(x: Double<f64>) -> Double<f64> {
            // GENERATE: sin_poly Double<f64> 9
            const K3: Double<f64> = Double::new(
                f64::from_bits(0xBFC5555555555555),
                f64::from_bits(0xBC65555552CA81AF),
            ); // -1.666666666666666666666666009500e-1
            const K5: Double<f64> = Double::new(
                f64::from_bits(0x3F81111111111111),
                f64::from_bits(0x3C0110E7CB9CA23D),
            ); // 8.333333333333333333329065980117e-3
            const K7: Double<f64> = Double::new(
                f64::from_bits(0xBF2A01A01A01A01A),
                f64::from_bits(0xBB572C66C6063195),
            ); // -1.984126984126984126029916072419e-4
            const K9: Double<f64> = Double::new(
                f64::from_bits(0x3EC71DE3A556C731),
                f64::from_bits(0x3B475F97DC1DE3C8),
            ); // 2.755731922398588019213703214996e-6
            const K11: Double<f64> = Double::new(
                f64::from_bits(0xBE5AE64567F53D3B),
                f64::from_bits(0xBAE00414E69D7B82),
            ); // -2.505210838543523223314376036023e-8
            const K13: Double<f64> = Double::new(
                f64::from_bits(0x3DE61246139A1624),
                f64::from_bits(0x3A819409726D3F63),
            ); // 1.605904383439241982126200763062e-10
            const K15: Double<f64> = Double::new(
                f64::from_bits(0xBD6AE7F3C61EDD83),
                f64::from_bits(0xB9FE33789227A734),
            ); // -7.647163171398073276859390202776e-13
            const K17: Double<f64> = Double::new(
                f64::from_bits(0x3CE9529972C79724),
                f64::from_bits(0x3990F4EA24C63034),
            ); // 2.811379344450285307153501050905e-15
            const K19: Double<f64> = Double::new(
                f64::from_bits(0xBC62D1434218BDB5),
                f64::from_bits(0xB8F94EABF8C41E64),
            ); // -8.160760491214478515701493028147e-18

            let x2 = x * x;
            let x3 = x2 * x;

            x + horner!(x3, x2, [K3, K5, K7, K9, K11, K13, K15, K17, K19])
        }

        fn cos(x: Double<f64>) -> Double<f64> {
            // GENERATE: cos_poly Double<f64> 8
            const K4: Double<f64> = Double::new(
                f64::from_bits(0x3FA5555555555555),
                f64::from_bits(0x3C4555424437931C),
            ); // 4.166666666666666666663512319225e-2
            const K6: Double<f64> = Double::new(
                f64::from_bits(0xBF56C16C16C16C16),
                f64::from_bits(0xBC0803A756A07D07),
            ); // -1.388888888888888887780931530771e-3
            const K8: Double<f64> = Double::new(
                f64::from_bits(0x3EFA01A01A01A015),
                f64::from_bits(0x3BA26CDCC55A9E8C),
            ); // 2.480158730158728657597925339703e-5
            const K10: Double<f64> = Double::new(
                f64::from_bits(0xBE927E4FB7789792),
                f64::from_bits(0xBB29D6F71F4B35F6),
            ); // -2.755731922397533319638598669759e-7
            const K12: Double<f64> = Double::new(
                f64::from_bits(0x3E21EED8EFE8FA9B),
                f64::from_bits(0x3ADE55AD93D3B7C7),
            ); // 2.087675698356730672524610350407e-9
            const K14: Double<f64> = Double::new(
                f64::from_bits(0xBDA9397481DD955F),
                f64::from_bits(0xBA566962E337915B),
            ); // -1.147074454371334149048113104850e-11
            const K16: Double<f64> = Double::new(
                f64::from_bits(0x3D2AE7BB51F3606D),
                f64::from_bits(0x39D134E8DEC3DCFC),
            ); // 4.779323963823227459005447965865e-14
            const K18: Double<f64> = Double::new(
                f64::from_bits(0xBCA6556752093B7B),
                f64::from_bits(0xB94DA7E65CAE91BA),
            ); // -1.549705349810074051809391126314e-16

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

            let x = Double::from(x);
            let k = Double::from(2.0);

            let d = x - k;
            let (i, i_f) = round_as_i_f(d.hi());
            let f = (d - i_f).normalize();

            // lgf = ln(Γ(k + f))
            let lgf = {
                // GENERATE: ln_gamma_poly Double<f64> 36 2 -0.5 0.50001
                const K1: Double<f64> = Double::new(
                    f64::from_bits(0x3FDB0EE6072093CE),
                    f64::from_bits(0x3C56CB90701FBFAB),
                ); // 4.227843350984671393934879099176e-1
                const K2: Double<f64> = Double::new(
                    f64::from_bits(0x3FD4A34CC4A60FA6),
                    f64::from_bits(0x3C71873D89122009),
                ); // 3.224670334241132182362075833230e-1
                const K3: Double<f64> = Double::new(
                    f64::from_bits(0xBFB13E001A557606),
                    f64::from_bits(0xBC6024BA0E83BFD5),
                ); // -6.735230105319809513324605383750e-2
                const K4: Double<f64> = Double::new(
                    f64::from_bits(0x3F951322AC7D8483),
                    f64::from_bits(0x3C3AFC89088D3161),
                ); // 2.058080842778454787900092414132e-2
                const K5: Double<f64> = Double::new(
                    f64::from_bits(0xBF7E404FC218F5F1),
                    f64::from_bits(0xBC20DACEC16EE7F0),
                ); // -7.385551028673985266273097139096e-3
                const K6: Double<f64> = Double::new(
                    f64::from_bits(0x3F67ADD6EADB6C2F),
                    f64::from_bits(0x3C1A921F5B04C168),
                ); // 2.890510330741523285752986798756e-3
                const K7: Double<f64> = Double::new(
                    f64::from_bits(0xBF538AC5C2BF8E07),
                    f64::from_bits(0xBC09D6CFCCB80A0B),
                ); // -1.192753911703260977113961817521e-3
                const K8: Double<f64> = Double::new(
                    f64::from_bits(0x3F40B36AF86396E8),
                    f64::from_bits(0x3BF7CB3CB056F839),
                ); // 5.096695247430424223358306313901e-4
                const K9: Double<f64> = Double::new(
                    f64::from_bits(0xBF2D3FD4C76D2FC7),
                    f64::from_bits(0xBBD8E033E72019D5),
                ); // -2.231547584535793797590581821663e-4
                const K10: Double<f64> = Double::new(
                    f64::from_bits(0x3F1A127B0F17D65A),
                    f64::from_bits(0x3BA9B6784172BC37),
                ); // 9.945751278180853370278043874578e-5
                const K11: Double<f64> = Double::new(
                    f64::from_bits(0xBF078DE5BD7C81EE),
                    f64::from_bits(0xBBBEFAB70E1142AB),
                ); // -4.492623673813314182985614684945e-5
                const K12: Double<f64> = Double::new(
                    f64::from_bits(0x3EF580DCEE66EB02),
                    f64::from_bits(0x3B92C4F616A00980),
                ); // 2.050721277567069206039175050237e-5
                const K13: Double<f64> = Double::new(
                    f64::from_bits(0xBEE3CBC963CE223F),
                    f64::from_bits(0xBB97BC59649CEB67),
                ); // -9.439488275268391195481996641389e-6
                const K14: Double<f64> = Double::new(
                    f64::from_bits(0x3ED2597A39F34A9A),
                    f64::from_bits(0x3B05EDEDC0CB9EC7),
                ); // 4.374866789907472930074701780064e-6
                const K15: Double<f64> = Double::new(
                    f64::from_bits(0xBEC11B2EB767965A),
                    f64::from_bits(0xBB7CC38C10C3A636),
                ); // -2.039215753801485578508317068122e-6
                const K16: Double<f64> = Double::new(
                    f64::from_bits(0x3EB0064CDEB234CD),
                    f64::from_bits(0x3B6311A0E35CF672),
                ); // 9.551412130410533460958461155424e-7
                const K17: Double<f64> = Double::new(
                    f64::from_bits(0xBE9E2600D93C5BBD),
                    f64::from_bits(0xBB42A9546E3E0B9F),
                ); // -4.492469198742686574396636298584e-7
                const K18: Double<f64> = Double::new(
                    f64::from_bits(0x3E8C76BBB3EDB719),
                    f64::from_bits(0x3B33DA9E9CA6101F),
                ); // 2.120718480507544680937399579938e-7
                const K19: Double<f64> = Double::new(
                    f64::from_bits(0xBE7AF5A6CBE1CB09),
                    f64::from_bits(0xBB169E3877737C16),
                ); // -1.004322482693896438213232902979e-7
                const K20: Double<f64> = Double::new(
                    f64::from_bits(0x3E699B93C28688CB),
                    f64::from_bits(0x3B139BDFC9EF1524),
                ); // 4.769810174893028071660751379732e-8
                const K21: Double<f64> = Double::new(
                    f64::from_bits(0xBE5862C72F692EEE),
                    f64::from_bits(0xBB0E316C2193A450),
                ); // -2.271109430578690504772161414684e-8
                const K22: Double<f64> = Double::new(
                    f64::from_bits(0x3E47469D9B61F70F),
                    f64::from_bits(0x3AE569F4185DAB32),
                ); // 1.083865873115697167381853909414e-8
                const K23: Double<f64> = Double::new(
                    f64::from_bits(0xBE36434B2CFAA05C),
                    f64::from_bits(0xBAE5DEA27A8DC6E3),
                ); // -5.183477383138976749176573867046e-9
                const K24: Double<f64> = Double::new(
                    f64::from_bits(0x3E2555AA48ABE538),
                    f64::from_bits(0x3AD5F6B7BDCFDBE3),
                ); // 2.483677768113273727948991063949e-9
                const K25: Double<f64> = Double::new(
                    f64::from_bits(0xBE147B07123F8A41),
                    f64::from_bits(0xBACF3F85A54829EA),
                ); // -1.192126460908977237512502955673e-9
                const K26: Double<f64> = Double::new(
                    f64::from_bits(0x3E03B138614E04CA),
                    f64::from_bits(0x3AA0DDC25F6E068E),
                ); // 5.731203867472705309374566811387e-10
                const K27: Double<f64> = Double::new(
                    f64::from_bits(0xBDF2F7A8BC265BEE),
                    f64::from_bits(0xBAAE7F2D7A4A4B56),
                ); // -2.760122652568227958083301686425e-10
                const K28: Double<f64> = Double::new(
                    f64::from_bits(0x3DE24B63D0396DD6),
                    f64::from_bits(0x3A972216AB936773),
                ); // 1.331099467734351145993975832208e-10
                const K29: Double<f64> = Double::new(
                    f64::from_bits(0xBDD19A226B1B342E),
                    f64::from_bits(0xBA818FC72D3B4DC2),
                ); // -6.403602195081112747904879736530e-11
                const K30: Double<f64> = Double::new(
                    f64::from_bits(0x3DC0F86EE2AA5CA7),
                    f64::from_bits(0x3A0EBBA7DC467B63),
                ); // 3.086905412985665043714690981573e-11
                const K31: Double<f64> = Double::new(
                    f64::from_bits(0xBDB10138281FCB26),
                    f64::from_bits(0xBA50DF44EB88D91E),
                ); // -1.546574197507270976808910881965e-11
                const K32: Double<f64> = Double::new(
                    f64::from_bits(0x3DA0C515DBA7298A),
                    f64::from_bits(0x3A51DBAF4725098D),
                ); // 7.626051582023145483915710453512e-12
                const K33: Double<f64> = Double::new(
                    f64::from_bits(0xBD890EE290F302EF),
                    f64::from_bits(0xBA118BEAC09B56F9),
                ); // -2.848781221698287440532080399828e-12
                const K34: Double<f64> = Double::new(
                    f64::from_bits(0x3D75F2B4F03A237B),
                    f64::from_bits(0x3A29857FE134A8FF),
                ); // 1.247603529502751232341366700243e-12
                const K35: Double<f64> = Double::new(
                    f64::from_bits(0xBD79330AEAD74119),
                    f64::from_bits(0xBA24A3F42B408743),
                ); // -1.432419215660462389905381982977e-12
                const K36: Double<f64> = Double::new(
                    f64::from_bits(0x3D6AAB78B8436D26),
                    f64::from_bits(0x39D09B400E90651D),
                ); // 7.580016127701884965269890955627e-13

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
                        v = v + 1.0;
                        gi = (gi * v).normalize();
                        if gi.hi().abs() > 1e150 {
                            gi = gi.pmul1(f64::exp2i_fast(-500));
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
                        gi = gi.pmul1(f64::exp2i_fast(500));
                        descale += 500;
                    }
                    for _ in 1..-i {
                        v = v + 1.0;
                        gi = (gi * v).normalize();
                        if gi.hi().abs() > 1e150 {
                            gi = gi.pmul1(f64::exp2i_fast(-500));
                            descale -= 500;
                        }
                    }
                    // ln(abs(Γ(x))) = ln(Γ(k + f)) - ln(abs(prod(j = 0 to |i| - 1, x + j))) = lgf - ln(abs(gi))
                    let sign = if gi.hi().sign() { -1 } else { 1 };
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
            let reflect = x.sign();
            let nx = if reflect {
                Double::new_sub11(1.0, x)
            } else {
                Double::from(x)
            };

            // p = P(1 / nx)
            let p = {
                // GENERATE: gamma_lanczos_poly Double<f64> 9 0.5 2.93e-39 0.022223
                const K0: Double<f64> = Double::new(
                    f64::from_bits(0x40040D931FF62705),
                    f64::from_bits(0x3CB2CAF948465E4E),
                ); // 2.506628274631000502415765307447e0
                const K1: Double<f64> = Double::new(
                    f64::from_bits(0x3FCABCC42A9D895C),
                    f64::from_bits(0x3C790E92BDA69D6F),
                ); // 2.088856895525833752011156633309e-1
                const K2: Double<f64> = Double::new(
                    f64::from_bits(0x3F81D32D71BE5B93),
                    f64::from_bits(0x3C269ED00C83FB3C),
                ); // 8.703570398024307585854645818948e-3
                const K3: Double<f64> = Double::new(
                    f64::from_bits(0xBF7B8792FC787D15),
                    f64::from_bits(0xBC1B63009BB23951),
                ); // -6.721090474030041162568244073157e-3
                const K4: Double<f64> = Double::new(
                    f64::from_bits(0xBF42D92340C82F31),
                    f64::from_bits(0xBBFA29E990B614B8),
                ); // -5.752012380652925138697500384312e-4
                const K5: Double<f64> = Double::new(
                    f64::from_bits(0x3F601985A6CD5224),
                    f64::from_bits(0x3BF4A8C31E99B8DE),
                ); // 1.965294874404459741276072510530e-3
                const K6: Double<f64> = Double::new(
                    f64::from_bits(0x3F26E8C0502F1D3A),
                    f64::from_bits(0x3BD68EF2C4708A60),
                ); // 1.747832066867157454773808239276e-4
                const K7: Double<f64> = Double::new(
                    f64::from_bits(0xBF5851F0FA6CDAE1),
                    f64::from_bits(0xBC0814F9BB47FB56),
                ); // -1.484380083366758742765523474202e-3
                const K8: Double<f64> = Double::new(
                    f64::from_bits(0xBF20D52DE91D4F7B),
                    f64::from_bits(0xBBA57C94D79C2952),
                ); // -1.284235518662599273329735872778e-4
                const K9: Double<f64> = Double::new(
                    f64::from_bits(0x3F611F449B9249DD),
                    f64::from_bits(0x3BC3637631223E74),
                ); // 2.090105056205644835111295019112e-3

                let inv_nx = nx.recip();
                K0 + horner!(inv_nx, inv_nx, [K1, K2, K3, K4, K5, K6, K7, K8, K9])
            };

            let ln_nx = ln(nx);

            // ln(Γ(nx)) = (nx - 0.5) * ln(nx) - nx + ln(P(1 / nx))
            //           = nx * (ln(nx) - 1) - 0.5 * ln(nx) + ln(P(1 / nx))
            // The second form is used to avoid overflow in intermediate operations.
            let lg_nx = (ln_nx - 1.0) * nx - ln_nx.pmul1(0.5) + ln(p);

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
                let sign = if sinpix.hi().sign() { -1 } else { 1 };

                (lgx, sign)
            } else {
                // ln(abs(Γ(x))) = ln(Γ(nx))
                (lg_nx, 1)
            }
        };

        if y.hi() == f64::INFINITY {
            (f64::INFINITY, sign)
        } else {
            (y.to_single(), sign)
        }
    }
}

fn reduce_half_revs(x: f64) -> (u8, f64) {
    let x_exp = x.exponent();
    if x_exp < -2 {
        // |x| < 0.25
        return (0, x);
    } else if x_exp > 52 + 2 {
        return (0, 0.0);
    }

    let x_abs = x.abs();
    let x_scale = x_abs * 2.0;
    let (xi, xi_f) = round_as_i_f(x_scale);
    let x_frac = x_scale - xi_f;

    let n = xi as u8;
    let y = x_frac * 0.5;

    if x.sign() {
        (n.wrapping_neg() & 3, -y)
    } else {
        (n & 3, y)
    }
}
