use super::f64x2::F64x2;
use crate::generic::{round_fi, scalbn};
use crate::traits::Float as _;

// GENERATE: consts f64 LOG2_E LOG2_10
const LOG2_E: f64 = f64::from_bits(0x3FF71547652B82FE); // 1.4426950408889634e0
const LOG2_10: f64 = f64::from_bits(0x400A934F0979A371); // 3.321928094887362e0

// GENERATE: consts F64x2 LN_2 LOG10_2 LN_10
const LN_2: F64x2 = F64x2::from_bits(0x3FE62E42FEFA39EF, 0x3C7ABC9E3B39803F); // 6.931471805599453094172321214582e-1
const LOG10_2: F64x2 = F64x2::from_bits(0x3FD34413509F79FF, 0xBC49DC1DA994FD21); // 3.010299956639811952137388947245e-1
const LN_10: F64x2 = F64x2::from_bits(0x40026BB1BBB55516, 0xBCAF48AD494EA3E9); // 2.302585092994045684017991454684e0

impl crate::generic::Exp for f64 {
    fn exp_finite(x: Self) -> Self {
        if x.exponent() >= 10 {
            if x.is_sign_negative() {
                return 0.0;
            } else {
                return f64::INFINITY;
            }
        }

        // x = k*ln(2) + r
        // k is an integer
        // |r| <= 0.5*ln(2)
        let (kf, k) = round_fi(x * LOG2_E);
        let r = x - kf * LN_2;

        // exp(r) = exp(r / s)^s
        // s = 2^5

        let rs = r.scalbn_fast(-5);

        // GENERATE: exp_m1_poly F64x2 6 -0.011 0.011
        const K2: F64x2 = F64x2::from_bits(0x3FE0000000000000, 0x3C5626A31E1B3BCF); // 5.000000000000000048032165218071e-1
        const K3: F64x2 = F64x2::from_bits(0x3FC5555555555555, 0x3C68E2B24638DCCF); // 1.666666666666666682071875965350e-1
        const K4: F64x2 = F64x2::from_bits(0x3FA555555554A281, 0xBC4822E2FBE19B24); // 4.166666666634899917392087166581e-2
        const K5: F64x2 = F64x2::from_bits(0x3F81111111108755, 0xBC13BCE0AFDAB175); // 8.333333333272166600425042119601e-3
        const K6: F64x2 = F64x2::from_bits(0x3F56C171BA3E471F, 0x3BD27D6517FDE5F9); // 1.388894140266612252425777542649e-3
        const K7: F64x2 = F64x2::from_bits(0x3F2A01A6678168B9, 0xBB8DBDB99C80076E); // 1.984134321471647986939491602879e-4

        let rs2 = rs.square();
        let mut t = rs + horner!(rs2, rs, [K2, K3, K4, K5, K6, K7]);

        // exp(r) - 1 = (t + 1)^s - 1
        // (t + 1)^2 - 1 = t^2 + 2*t
        for _ in 0..5 {
            t = t.square() + t.twice();
        }

        // y = exp(x) = exp(r) * 2^k = (t + 1) * 2^k
        (t + 1.0).scalbn_to_f64(k as i32)
    }

    fn exp_m1_finite(x: Self) -> Self {
        let x_exp = x.exponent();
        if x_exp >= 10 {
            if x.is_sign_negative() {
                return -1.0;
            } else {
                return f64::INFINITY;
            }
        } else if x_exp <= -600 {
            return x;
        }

        // x = k*ln(2) + r
        // k is an integer
        // |r| <= 0.5*ln(2)
        let (kf, k) = round_fi(x * LOG2_E);
        if k <= -55 {
            return -1.0;
        }
        let k = k as i32;
        let r = x - kf * LN_2;

        // exp(r) = exp(r / s)^s
        // s = 2^5

        let rs = r.scalbn_fast(-5);

        // GENERATE: exp_m1_poly F64x2 6 -0.011 0.011
        const K2: F64x2 = F64x2::from_bits(0x3FE0000000000000, 0x3C5626A31E1B3BCF); // 5.000000000000000048032165218071e-1
        const K3: F64x2 = F64x2::from_bits(0x3FC5555555555555, 0x3C68E2B24638DCCF); // 1.666666666666666682071875965350e-1
        const K4: F64x2 = F64x2::from_bits(0x3FA555555554A281, 0xBC4822E2FBE19B24); // 4.166666666634899917392087166581e-2
        const K5: F64x2 = F64x2::from_bits(0x3F81111111108755, 0xBC13BCE0AFDAB175); // 8.333333333272166600425042119601e-3
        const K6: F64x2 = F64x2::from_bits(0x3F56C171BA3E471F, 0x3BD27D6517FDE5F9); // 1.388894140266612252425777542649e-3
        const K7: F64x2 = F64x2::from_bits(0x3F2A01A6678168B9, 0xBB8DBDB99C80076E); // 1.984134321471647986939491602879e-4

        let rs2 = rs.square();
        let mut t = rs + horner!(rs2, rs, [K2, K3, K4, K5, K6, K7]);

        // exp(r) - 1 = (t + 1)^s - 1
        // (t + 1)^2 - 1 = t^2 + 2*t
        for _ in 0..5 {
            t = t.square() + t.twice();
        }

        // y = exp(x) - 1 = exp(r) * 2^k - 1 = 2^k * (t + (1 - 2^-k))
        (t + F64x2::sub11(1.0, scalbn(1.0, -k))).scalbn_to_f64(k)
    }

    fn exp2_finite(x: Self) -> Self {
        if x.exponent() >= 11 {
            if x.is_sign_negative() {
                return 0.0;
            } else {
                return f64::INFINITY;
            }
        }

        // x = k + r*log2(e)
        // k is an integer
        // |r| <= 0.5*ln(2)
        let (kf, k) = round_fi(x);
        let r = (x - kf) * LN_2;

        // exp(r) = exp(r / s)^s
        // s = 2^5

        let rs = r.scalbn_fast(-5);

        // GENERATE: exp_m1_poly F64x2 6 -0.011 0.011
        const K2: F64x2 = F64x2::from_bits(0x3FE0000000000000, 0x3C5626A31E1B3BCF); // 5.000000000000000048032165218071e-1
        const K3: F64x2 = F64x2::from_bits(0x3FC5555555555555, 0x3C68E2B24638DCCF); // 1.666666666666666682071875965350e-1
        const K4: F64x2 = F64x2::from_bits(0x3FA555555554A281, 0xBC4822E2FBE19B24); // 4.166666666634899917392087166581e-2
        const K5: F64x2 = F64x2::from_bits(0x3F81111111108755, 0xBC13BCE0AFDAB175); // 8.333333333272166600425042119601e-3
        const K6: F64x2 = F64x2::from_bits(0x3F56C171BA3E471F, 0x3BD27D6517FDE5F9); // 1.388894140266612252425777542649e-3
        const K7: F64x2 = F64x2::from_bits(0x3F2A01A6678168B9, 0xBB8DBDB99C80076E); // 1.984134321471647986939491602879e-4

        let rs2 = rs.square();
        let mut t = rs + horner!(rs2, rs, [K2, K3, K4, K5, K6, K7]);

        // exp(r) - 1 = (t + 1)^s - 1
        // (t + 1)^2 - 1 = t^2 + 2*t
        for _ in 0..5 {
            t = t.square() + t.twice();
        }

        // y = exp(x) = exp(r) * 2^k = (t + 1) * 2^k
        (t + 1.0).scalbn_to_f64(k as i32)
    }

    fn exp10_finite(x: Self) -> Self {
        if x.exponent() >= 9 {
            if x.is_sign_negative() {
                return 0.0;
            } else {
                return f64::INFINITY;
            }
        }

        // x = k*log10(2) + r*log10(e)
        // k is an integer
        // |r| <= 0.5*ln(2)
        let (kf, k) = round_fi(x * LOG2_10);
        let r = (x - kf * LOG10_2) * LN_10;

        // exp(r) = exp(r / s)^s
        // s = 2^5

        let rs = r.scalbn_fast(-5);

        // GENERATE: exp_m1_poly F64x2 6 -0.011 0.011
        const K2: F64x2 = F64x2::from_bits(0x3FE0000000000000, 0x3C5626A31E1B3BCF); // 5.000000000000000048032165218071e-1
        const K3: F64x2 = F64x2::from_bits(0x3FC5555555555555, 0x3C68E2B24638DCCF); // 1.666666666666666682071875965350e-1
        const K4: F64x2 = F64x2::from_bits(0x3FA555555554A281, 0xBC4822E2FBE19B24); // 4.166666666634899917392087166581e-2
        const K5: F64x2 = F64x2::from_bits(0x3F81111111108755, 0xBC13BCE0AFDAB175); // 8.333333333272166600425042119601e-3
        const K6: F64x2 = F64x2::from_bits(0x3F56C171BA3E471F, 0x3BD27D6517FDE5F9); // 1.388894140266612252425777542649e-3
        const K7: F64x2 = F64x2::from_bits(0x3F2A01A6678168B9, 0xBB8DBDB99C80076E); // 1.984134321471647986939491602879e-4

        let rs2 = rs.square();
        let mut t = rs + horner!(rs2, rs, [K2, K3, K4, K5, K6, K7]);

        // exp(r) - 1 = (t + 1)^s - 1
        // (t + 1)^2 - 1 = t^2 + 2*t
        for _ in 0..5 {
            t = t.square() + t.twice();
        }

        // y = exp(x) = exp(r) * 2^k = (t + 1) * 2^k
        (t + 1.0).scalbn_to_f64(k as i32)
    }
}
