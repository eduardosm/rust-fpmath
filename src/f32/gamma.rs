use super::log_core::log_core;
use crate::double::Double;
use crate::generic::round_as_i_f;
use crate::scalbn;
use crate::traits::Float as _;

// GENERATE: consts Double<f32> PI LN_2 LN_PI
const PI: Double<f32> = Double::new(f32::from_bits(0x40490FDA), f32::from_bits(0x34222169)); // 3.1415926535898e0
const LN_2: Double<f32> = Double::new(f32::from_bits(0x3F317217), f32::from_bits(0x3377D1CF)); // 6.9314718055995e-1
const LN_PI: Double<f32> = Double::new(f32::from_bits(0x3F928682), f32::from_bits(0x330E7A1C)); // 1.1447298858494e0

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

        let x = Double::from(x);
        let k = Double::from(2.875);

        let d = x - k;
        let (i, i_f) = round_as_i_f(d.hi());
        let f = (d - i_f).normalize();

        // gf = Γ(k + f)
        let gf = {
            // GENERATE: gamma_poly Double<f32> 16 2.875 -0.5 0.5
            const K0: Double<f32> =
                Double::new(f32::from_bits(0x3FE4D3B5), f32::from_bits(0x33F2BB89)); // 1.7877108988969e0
            const K1: Double<f32> =
                Double::new(f32::from_bits(0x3FC793AA), f32::from_bits(0x335DCF91)); // 1.5591939012080e0
            const K2: Double<f32> =
                Double::new(f32::from_bits(0x3F8688C8), f32::from_bits(0x33CA4A0C)); // 1.0510493266812e0
            const K3: Double<f32> =
                Double::new(f32::from_bits(0x3EF0FA16), f32::from_bits(0x3270AE61)); // 4.7065801829339e-1
            const K4: Double<f32> =
                Double::new(f32::from_bits(0x3E4159AC), f32::from_bits(0x31A3AE90)); // 1.8881863832014e-1
            const K5: Double<f32> =
                Double::new(f32::from_bits(0x3D70F959), f32::from_bits(0x31065CD1)); // 5.8831548410854e-2
            const K6: Double<f32> =
                Double::new(f32::from_bits(0x3C9207B6), f32::from_bits(0x309E98EE)); // 1.7825943640472e-2
            const K7: Double<f32> =
                Double::new(f32::from_bits(0x3B8A9165), f32::from_bits(0x2CD502E2)); // 4.2287581660790e-3
            const K8: Double<f32> =
                Double::new(f32::from_bits(0x3A8FE804), f32::from_bits(0x2DB5C573)); // 1.0979180430890e-3
            const K9: Double<f32> =
                Double::new(f32::from_bits(0x394C0448), f32::from_bits(0x2D415F0B)); // 1.9456552043785e-4
            const K10: Double<f32> =
                Double::new(f32::from_bits(0x3859FA03), f32::from_bits(0x2C34A870)); // 5.1969675789933e-5
            const K11: Double<f32> =
                Double::new(f32::from_bits(0x36A4EBD1), f32::from_bits(0x2AE360B1)); // 4.9150339183116e-6
            const K12: Double<f32> =
                Double::new(f32::from_bits(0x36241572), f32::from_bits(0x29EA5DA6)); // 2.4450388213485e-6
            const K13: Double<f32> =
                Double::new(f32::from_bits(0xB41DD711), f32::from_bits(0xA83449F9)); // -1.4700006125073e-7
            const K14: Double<f32> =
                Double::new(f32::from_bits(0x342F0A78), f32::from_bits(0x286B1839)); // 1.6301954869872e-7
            const K15: Double<f32> =
                Double::new(f32::from_bits(0xB348562D), f32::from_bits(0xA7154B6F)); // -4.6644507219122e-8
            const K16: Double<f32> =
                Double::new(f32::from_bits(0x32A9DBE6), f32::from_bits(0x26E45280)); // 1.9774189204973e-8

            K0 + horner!(
                f,
                f,
                [
                    K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16
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
                    if gi.hi().abs() > 1e20 {
                        gi = gi.pmul1(f32::exp2i_fast(-50));
                        descale += 50;
                    }
                }
                // Γ(x) = Γ(k + f) * prod(j = 0 to i - 1, k + f + j) = gf * gi
                gf * gi
            }
            core::cmp::Ordering::Less => {
                // gi = prod(j = 0 to |i| - 1, x + j)
                let mut v = x;
                let mut gi = v;
                if gi.hi().abs() < 1e-20 {
                    gi = gi.pmul1(f32::exp2i_fast(50));
                    descale += 50;
                }
                for _ in 1..-i {
                    v = v + 1.0;
                    gi = (gi * v).normalize();
                    if gi.hi().abs() > 1e20 {
                        gi = gi.pmul1(f32::exp2i_fast(-50));
                        descale -= 50;
                    }
                }
                // Γ(x) = Γ(k + f) / prod(j = 0 to |i| - 1, x + j) = gf / gi
                gf / gi
            }
        };

        scalbn(y.to_single(), descale)
    }

    fn ln_gamma_finite(x: Self) -> (Self, i8) {
        fn ln(x: Double<f32>) -> Double<f32> {
            let xm1 = (x - 1.0).normalize();
            if xm1.hi().abs() < 0.0002 {
                // GENERATE: ln_1p_poly Double<f32> 4 -0.00021 0.00021
                const K2: Double<f32> =
                    Double::new(f32::from_bits(0xBEFFFFFF), f32::from_bits(0xB3000000)); // -5.0000000000000e-1
                const K3: Double<f32> =
                    Double::new(f32::from_bits(0x3EAAAAAA), f32::from_bits(0x32AAAAAB)); // 3.3333333333333e-1
                const K4: Double<f32> =
                    Double::new(f32::from_bits(0xBE800000), f32::from_bits(0xB21DD65B)); // -2.5000000918734e-1
                const K5: Double<f32> =
                    Double::new(f32::from_bits(0x3E4CCCCD), f32::from_bits(0x31FFBE70)); // 2.0000001042336e-1

                let xm1_2 = xm1.square();
                xm1 + horner!(xm1_2, xm1, [K2, K3, K4, K5])
            } else {
                // GENERATE: ln_1p_poly Double<f32> 7 -0.0001 0.0625
                const K2: Double<f32> =
                    Double::new(f32::from_bits(0xBEFFFFFF), f32::from_bits(0xB2FFFE14)); // -4.9999999999913e-1
                const K3: Double<f32> =
                    Double::new(f32::from_bits(0x3EAAAAAA), f32::from_bits(0x32A794DA)); // 3.3333333297417e-1
                const K4: Double<f32> =
                    Double::new(f32::from_bits(0xBE7FFFFC), f32::from_bits(0xB224B7B6)); // -2.4999994998318e-1
                const K5: Double<f32> =
                    Double::new(f32::from_bits(0x3E4CCBEE), f32::from_bits(0x3235DD62)); // 1.9999669060721e-1
                const K6: Double<f32> =
                    Double::new(f32::from_bits(0xBE2A8BD6), f32::from_bits(0xB2706478)); // -1.6654907076110e-1
                const K7: Double<f32> =
                    Double::new(f32::from_bits(0x3E0FE9E5), f32::from_bits(0x30F1D782)); // 1.4054067608843e-1
                const K8: Double<f32> =
                    Double::new(f32::from_bits(0xBDCE7E72), f32::from_bits(0xB1BEC5BD)); // -1.0082711834613e-1

                let (k, lo, ln_hi) = log_core(x);
                let lo2 = lo * lo;
                let ln_lo = lo + horner!(lo2, lo, [K2, K3, K4, K5, K6, K7, K8]);

                // ln(x) = ln(2^k * m) = k * ln(2) + ln(m)
                LN_2 * k + (ln_hi + ln_lo)
            }
        }

        fn sin(x: Double<f32>) -> Double<f32> {
            // GENERATE: sin_poly Double<f32> 5
            const K3: Double<f32> =
                Double::new(f32::from_bits(0xBE2AAAAA), f32::from_bits(0xB22AA91D)); // -1.6666666666631e-1
            const K5: Double<f32> =
                Double::new(f32::from_bits(0x3C088888), f32::from_bits(0x30064E28)); // 8.3333333252278e-3
            const K7: Double<f32> =
                Double::new(f32::from_bits(0xB9500CFC), f32::from_bits(0xAD28F32F)); // -1.9841263798241e-4
            const K9: Double<f32> =
                Double::new(f32::from_bits(0x3638EBBC), f32::from_bits(0x2909A447)); // 2.7555352367962e-6
            const K11: Double<f32> =
                Double::new(f32::from_bits(0xB2D4B292), f32::from_bits(0xA6273E8A)); // -2.4761253184112e-8

            let x2 = x * x;
            let x3 = x2 * x;

            x + horner!(x3, x2, [K3, K5, K7, K9, K11])
        }

        fn cos(x: Double<f32>) -> Double<f32> {
            // GENERATE: cos_poly Double<f32> 5
            const K4: Double<f32> =
                Double::new(f32::from_bits(0x3D2AAAAA), f32::from_bits(0x312AA986)); // 4.1666666666602e-2
            const K6: Double<f32> =
                Double::new(f32::from_bits(0xBAB60B60), f32::from_bits(0xAEB3B22B)); // -1.3888888878209e-3
            const K8: Double<f32> =
                Double::new(f32::from_bits(0x37D00CFD), f32::from_bits(0x2B2234A4)); // 2.4801580942599e-5
            const K10: Double<f32> =
                Double::new(f32::from_bits(0xB493F013), f32::from_bits(0xA845AAA2)); // -2.7555561774953e-7
            const K12: Double<f32> =
                Double::new(f32::from_bits(0x310DE345), f32::from_bits(0x244EFC36)); // 2.0647388692416e-9

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

            let x = Double::from(x);
            let k = Double::from(2.0);

            let d = x - k;
            let (i, i_f) = round_as_i_f(d.hi());
            let f = (d - i_f).normalize();

            // lgf = ln(Γ(k + f))
            let lgf = {
                // GENERATE: ln_gamma_poly Double<f32> 15 2 -0.5 0.50001
                const K1: Double<f32> =
                    Double::new(f32::from_bits(0x3ED87730), f32::from_bits(0x31E4127D)); // 4.2278433509847e-1
                const K2: Double<f32> =
                    Double::new(f32::from_bits(0x3EA51A66), f32::from_bits(0x3194C281)); // 3.2246703342418e-1
                const K3: Double<f32> =
                    Double::new(f32::from_bits(0xBD89F000), f32::from_bits(0xB1D2B147)); // -6.7352301053834e-2
                const K4: Double<f32> =
                    Double::new(f32::from_bits(0x3CA89915), f32::from_bits(0x30453FA9)); // 2.0580808418340e-2
                const K5: Double<f32> =
                    Double::new(f32::from_bits(0xBBF2027D), f32::from_bits(0xAFF8F494)); // -7.3855509853372e-3
                const K6: Double<f32> =
                    Double::new(f32::from_bits(0x3B3D6EB9), f32::from_bits(0x2DD5EECB)); // 2.8905107417282e-3
                const K7: Double<f32> =
                    Double::new(f32::from_bits(0xBA9C5637), f32::from_bits(0xAEC7DFCA)); // -1.1927550403353e-3
                const K8: Double<f32> =
                    Double::new(f32::from_bits(0x3A059AD0), f32::from_bits(0x2E7E70E8)); // 5.0966168019864e-4
                const K9: Double<f32> =
                    Double::new(f32::from_bits(0xB969FAC8), f32::from_bits(0xAD63DFCE)); // -2.2314036163347e-4
                const K10: Double<f32> =
                    Double::new(f32::from_bits(0x38D0BCDB), f32::from_bits(0x2C16FF19)); // 9.9533901774812e-5
                const K11: Double<f32> =
                    Double::new(f32::from_bits(0xB83CD94D), f32::from_bits(0xAC7AE746)); // -4.5025073559705e-5
                const K12: Double<f32> =
                    Double::new(f32::from_bits(0x37A8B7C3), f32::from_bits(0x2B50EE03)); // 2.0112732105310e-5
                const K13: Double<f32> =
                    Double::new(f32::from_bits(0xB7182E46), f32::from_bits(0xAA8C536F)); // -9.0706801293655e-6
                const K14: Double<f32> =
                    Double::new(f32::from_bits(0x36B4D890), f32::from_bits(0x2A6E5B20)); // 5.3896294346209e-6
                const K15: Double<f32> =
                    Double::new(f32::from_bits(0xB636F594), f32::from_bits(0xAA1B0EEE)); // -2.7263060032275e-6

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
                    let mut descale = 0i16;
                    let mut v = k + f;
                    let mut gi = v;
                    for _ in 1..i {
                        v = v + 1.0;
                        gi = (gi * v).normalize();
                        if gi.hi().abs() > 1e20 {
                            gi = gi.pmul1(f32::exp2i_fast(-50));
                            descale += 50;
                        }
                    }
                    // ln(abs(Γ(x))) = ln(Γ(k + f)) + ln(prod(j = 0 to i - 1, k + f + j)) = lgf + ln(gi)
                    (lgf + ln(gi) + LN_2 * f32::from(descale), 1)
                }
                core::cmp::Ordering::Less => {
                    // gi = prod(j = 0 to |i| - 1, x + j)
                    let mut descale = 0i16;
                    let mut v = x;
                    let mut gi = v;
                    if gi.hi().abs() < 1e-20 {
                        gi = gi.pmul1(f32::exp2i_fast(50));
                        descale += 50;
                    }
                    for _ in 1..-i {
                        v = v + 1.0;
                        gi = (gi * v).normalize();
                        if gi.hi().abs() > 1e20 {
                            gi = gi.pmul1(f32::exp2i_fast(-50));
                            descale -= 50;
                        }
                    }
                    // ln(abs(Γ(x))) = ln(Γ(k + f)) - ln(abs(prod(j = 0 to |i| - 1, x + j))) = lgf - ln(abs(gi))
                    let sign = if gi.hi().sign() { -1 } else { 1 };
                    (lgf - ln(gi.abs()) + LN_2 * f32::from(descale), sign)
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
                // GENERATE: gamma_lanczos_poly Double<f32> 4 0.5 2.93e-39 0.022223
                const K0: Double<f32> =
                    Double::new(f32::from_bits(0x40206C98), f32::from_bits(0x347FB13A)); // 2.5066282746310e0
                const K1: Double<f32> =
                    Double::new(f32::from_bits(0x3E55E621), f32::from_bits(0x31A83B58)); // 2.0888568950560e-1
                const K2: Double<f32> =
                    Double::new(f32::from_bits(0x3C0E997D), f32::from_bits(0x3035D007)); // 8.7035873068547e-3
                const K3: Double<f32> =
                    Double::new(f32::from_bits(0xBBDC4E75), f32::from_bits(0xAF83C2AF)); // -6.7232202341265e-3
                const K4: Double<f32> =
                    Double::new(f32::from_bits(0xB9F42FCE), f32::from_bits(0xADC8EDDC)); // -4.6574926192638e-4

                let inv_nx = nx.recip();
                K0 + horner!(inv_nx, inv_nx, [K1, K2, K3, K4])
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

        if y.hi() == f32::INFINITY {
            (f32::INFINITY, sign)
        } else {
            (y.to_single(), sign)
        }
    }
}

fn reduce_half_revs(x: f32) -> (u8, f32) {
    let x_exp = x.exponent();
    if x_exp < -2 {
        // |x| < 0.25
        return (0, x);
    } else if x_exp > 23 + 2 {
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
