use super::log_core::log_core_f64;
use crate::generic::{reduce_half_revs, round_fi};
use crate::traits::Float as _;

// GENERATE: consts f64 PI LN_2 LN_PI
const PI: f64 = f64::from_bits(0x400921FB54442D18); // 3.141592653589793e0
const LN_2: f64 = f64::from_bits(0x3FE62E42FEFA39EF); // 6.931471805599453e-1
const LN_PI: f64 = f64::from_bits(0x3FF250D048E7A1BD); // 1.1447298858494002e0

impl crate::generic::Gamma for f32 {
    fn gamma_finite(x: Self) -> Self {
        if x >= 36.0 {
            // Overflow
            return f32::INFINITY;
        } else if x <= -43.0 {
            // Underflow
            let xi = (-x) as u32;
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

        let x = f64::from(x);
        let k = 2.875;

        let d = x - k;
        let (i_f, i) = round_fi(d);
        let f = d - i_f;

        // gf = Γ(k + f)
        let gf = {
            // GENERATE: gamma_poly f64 10 2.875 -0.5 0.5
            const K0: f64 = f64::from_bits(0x3FFC9A76BE5776DF); // 1.7877108988972663e0
            const K1: f64 = f64::from_bits(0x3FF8F2754DE034AA); // 1.5591939012549951e0
            const K2: f64 = f64::from_bits(0x3FF0D11919467E31); // 1.0510493266409748e0
            const K3: f64 = f64::from_bits(0x3FDE1F42CAE80255); // 4.70658014441175e-1
            const K4: f64 = f64::from_bits(0x3FC82B358A50087A); // 1.8881863835855822e-1
            const K5: f64 = f64::from_bits(0x3FAE1F2E230C7F02); // 5.883163621751743e-2
            const K6: f64 = f64::from_bits(0x3F9240F86283980F); // 1.782596684724785e-2
            const K7: f64 = f64::from_bits(0x3F715151FAE5937A); // 4.227943644171778e-3
            const K8: f64 = f64::from_bits(0x3F51FBA622CBD87B); // 1.0975954457399784e-3
            const K9: f64 = f64::from_bits(0x3F29EF7115EF3777); // 1.9787078323790414e-4
            const K10: f64 = f64::from_bits(0x3F0C0EEDBA051121); // 5.3516988602927156e-5

            K0 + horner!(f, f, [K1, K2, K3, K4, K5, K6, K7, K8, K9, K10])
        };

        let y = match i.cmp(&0) {
            core::cmp::Ordering::Equal => gf,
            core::cmp::Ordering::Greater => {
                // gi = prod(j = 0 to i - 1, k + f + j)
                let mut v = k + f;
                let mut gi = v;
                for _ in 1..i {
                    v += 1.0;
                    gi *= v;
                }
                // Γ(x) = Γ(k + f) * prod(j = 0 to i - 1, k + f + j) = gf * gi
                gf * gi
            }
            core::cmp::Ordering::Less => {
                // gi = prod(j = 0 to |i| - 1, x + j)
                let mut v = x;
                let mut gi = v;
                for _ in 1..-i {
                    v += 1.0;
                    gi *= v;
                }
                // Γ(x) = Γ(k + f) / prod(j = 0 to |i| - 1, x + j) = gf / gi
                gf / gi
            }
        };

        y as f32
    }

    fn ln_gamma_finite(x: Self) -> (Self, i8) {
        fn ln(x: f64) -> f64 {
            let xm1 = x - 1.0;
            if xm1.exponent() <= -10 {
                // GENERATE: ln_1p_poly f64 4 -0.00196 0.00196
                const K2: f64 = f64::from_bits(0xBFDFFFFFFFFFC9F1); // -4.999999999992318e-1
                const K3: f64 = f64::from_bits(0x3FD555555554F43F); // 3.3333333333195364e-1
                const K4: f64 = f64::from_bits(0xBFD000035B35DA72); // -2.5000080020200877e-1
                const K5: f64 = f64::from_bits(0x3FC999A1376977B1); // 2.0000090795195782e-1

                let xm1_2 = xm1 * xm1;
                xm1 + horner!(xm1_2, xm1, [K2, K3, K4, K5])
            } else {
                // GENERATE: ln_1p_poly f64 6 -0.032 0.032
                const K2: f64 = f64::from_bits(0xBFE000000001FF45); // -5.000000000145312e-1
                const K3: f64 = f64::from_bits(0x3FD55555555FA9F6); // 3.3333333337091575e-1
                const K4: f64 = f64::from_bits(0xBFCFFFFF0A7CC9E0); // -2.4999988567431242e-1
                const K5: f64 = f64::from_bits(0x3FC999981E53B5B4); // 1.9999982338724254e-1
                const K6: f64 = f64::from_bits(0xBFC55CAB5870E137); // -1.6689054315953353e-1
                const K7: f64 = f64::from_bits(0x3FC25159F6260532); // 1.431076480767302e-1

                let (k, lo, ln_hi) = log_core_f64(x, 0);
                let lo2 = lo * lo;
                let ln_lo = lo + horner!(lo2, lo, [K2, K3, K4, K5, K6, K7]);

                // ln(x) = ln(2^k * m) = k * ln(2) + ln(m)
                LN_2 * k + (ln_hi + ln_lo)
            }
        }

        fn sin(x: f64) -> f64 {
            // GENERATE: sin_poly f64 5
            const K3: f64 = f64::from_bits(0xBFC55555555523A3); // -1.6666666666631355e-1
            const K5: f64 = f64::from_bits(0x3F81111110C9C505); // 8.333333325227796e-3
            const K7: f64 = f64::from_bits(0xBF2A019F951E65D9); // -1.9841263798241088e-4
            const K9: f64 = f64::from_bits(0x3EC71D77844D2238); // 2.7555352367962157e-6
            const K11: f64 = f64::from_bits(0xBE5A96524A73E8A7); // -2.476125318411175e-8

            let x2 = x * x;
            let x3 = x2 * x;

            x + horner!(x3, x2, [K3, K5, K7, K9, K11])
        }

        fn cos(x: f64) -> f64 {
            // GENERATE: cos_poly f64 5
            const K4: f64 = f64::from_bits(0x3FA55555555530CC); // 4.1666666666601765e-2
            const K6: f64 = f64::from_bits(0xBF56C16C1676455D); // -1.388888887820925e-3
            const K8: f64 = f64::from_bits(0x3EFA019FAA234A38); // 2.4801580942599414e-5
            const K10: f64 = f64::from_bits(0xBE927E026C5AAA26); // -2.755556177495311e-7
            const K12: f64 = f64::from_bits(0x3E21BC68A677E1B1); // 2.0647388692416302e-9

            let x2 = x * x;
            let x4 = x2 * x2;

            let t = horner!(x4, x2, [K4, K6, K8, K10, K12]);
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

            let x = f64::from(x);
            let k = 2.0;

            let d = x - k;
            let (i_f, i) = round_fi(d);
            let f = d - i_f;

            // lgf = ln(Γ(k + f))
            let lgf = {
                // GENERATE: ln_gamma_poly f64 15 2 -0.5 0.50001
                const K1: f64 = f64::from_bits(0x3FDB0EE6072093EA); // 4.227843350984687e-1
                const K2: f64 = f64::from_bits(0x3FD4A34CC4A6140B); // 3.2246703342417565e-1
                const K3: f64 = f64::from_bits(0xBFB13E001A5628EC); // -6.735230105383366e-2
                const K4: f64 = f64::from_bits(0x3F951322AC53FA8A); // 2.058080841833968e-2
                const K5: f64 = f64::from_bits(0xBF7E404FBF1E9289); // -7.385550985337227e-3
                const K6: f64 = f64::from_bits(0x3F67ADD72357BB2D); // 2.890510741728211e-3
                const K7: f64 = f64::from_bits(0xBF538AC6F8FBF93C); // -1.1927550403352935e-3
                const K8: f64 = f64::from_bits(0x3F40B35A1FCE1D04); // 5.096616801986424e-4
                const K9: f64 = f64::from_bits(0xBF2D3F591C7BF9B1); // -2.231403616334677e-4
                const K10: f64 = f64::from_bits(0x3F1A179B696FF18D); // 9.953390177481186e-5
                const K11: f64 = f64::from_bits(0xBF079B29BF5CE8B6); // -4.502507355970532e-5
                const K12: f64 = f64::from_bits(0x3EF516F86D0EE035); // 2.0112732105309796e-5
                const K13: f64 = f64::from_bits(0xBEE305C8C8C536EB); // -9.070680129365463e-6
                const K14: f64 = f64::from_bits(0x3ED69B120EE5B1FE); // 5.389629434620917e-6
                const K15: f64 = f64::from_bits(0xBEC6DEB29361DDB0); // -2.7263060032274846e-6

                horner!(
                    f,
                    f,
                    [
                        K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15
                    ]
                )
            };

            match i.cmp(&0) {
                core::cmp::Ordering::Equal => (lgf, 1),
                core::cmp::Ordering::Greater => {
                    // gi = prod(j = 0 to i - 1, k + f + j)
                    let mut v = k + f;
                    let mut gi = v;
                    for _ in 1..i {
                        v += 1.0;
                        gi *= v;
                    }
                    // ln(abs(Γ(x))) = ln(Γ(k + f)) + ln(prod(j = 0 to i - 1, k + f + j)) = lgf + ln(gi)
                    (lgf + ln(gi), 1)
                }
                core::cmp::Ordering::Less => {
                    // gi = prod(j = 0 to |i| - 1, x + j)
                    let mut v = x;
                    let mut gi = v;
                    for _ in 1..-i {
                        v += 1.0;
                        gi *= v;
                    }
                    // ln(abs(Γ(x))) = ln(Γ(k + f)) - ln(abs(prod(j = 0 to |i| - 1, x + j))) = lgf - ln(abs(gi))
                    let sign = if gi.is_sign_negative() { -1 } else { 1 };
                    (lgf - ln(gi.abs()), sign)
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
                1.0 - f64::from(x)
            } else {
                f64::from(x)
            };

            // p = P(1 / nx)
            let p = {
                // GENERATE: gamma_lanczos_poly f64 4 0.5 2.93e-39 0.022223
                const K0: f64 = f64::from_bits(0x40040D931FF62735); // 2.5066282746310216e0
                const K1: f64 = f64::from_bits(0x3FCABCC42A83B582); // 2.0888568950560332e-1
                const K2: f64 = f64::from_bits(0x3F81D32FB6BA00D0); // 8.703587306854749e-3
                const K3: f64 = f64::from_bits(0xBF7B89CEB07855D3); // -6.723220234126487e-3
                const K4: f64 = f64::from_bits(0xBF3E85F9D91DBB89); // -4.657492619263783e-4

                let inv_nx = nx.recip();
                K0 + horner!(inv_nx, inv_nx, [K1, K2, K3, K4])
            };

            // ln(Γ(nx)) = (nx - 0.5) * ln(nx) - nx + ln(P(1 / nx))
            let lg_nx = (nx - 0.5) * ln(nx) - nx + ln(p);

            if reflect {
                let (n, z) = reduce_half_revs(x);
                let z_rad = PI * f64::from(z);
                let sinpix = match n {
                    0 => sin(z_rad),
                    1 => cos(z_rad),
                    2 => -sin(z_rad),
                    3 => -cos(z_rad),
                    _ => unreachable!(),
                };

                // ln(abs(Γ(x))) = ln(π) - ln(abs(sin(πx))) - ln(Γ(1-x))
                let lgx = LN_PI - ln(sinpix.abs()) - lg_nx;
                let sign = if sinpix.is_sign_negative() { -1 } else { 1 };

                (lgx, sign)
            } else {
                // ln(abs(Γ(x))) = ln(Γ(nx))
                (lg_nx, 1)
            }
        };

        (y as f32, sign)
    }
}
