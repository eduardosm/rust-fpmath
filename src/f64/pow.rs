use super::f64x2::F64x2;
use super::log_core::log_core_f64;
use crate::generic::round_fi;
use crate::traits::Float as _;

// GENERATE: consts f64 LOG2_E
const LOG2_E: f64 = f64::from_bits(0x3FF71547652B82FE); // 1.4426950408889634e0

// GENERATE: consts F64x2 LN_2
const LN_2: F64x2 = F64x2::from_bits(0x3FE62E42FEFA39EF, 0x3C7ABC9E3B39803F); // 6.931471805599453094172321214582e-1

impl crate::generic::Pow for f64 {
    fn pow_finite(x: Self, xedelta: Self::Exp, y: Self, sign: bool) -> Self {
        pow_core(x, xedelta, y, sign)
    }

    fn powi_finite(x: Self, xedelta: Self::Exp, y: i32) -> Self {
        pow_core(x, xedelta, y.into(), x.is_sign_negative() && (y & 1) != 0)
    }
}

fn pow_core(x: f64, xedelta: i16, y: f64, sign: bool) -> f64 {
    #[inline]
    fn ln(x: f64, edelta: i16) -> F64x2 {
        // GENERATE: ln_1p_poly F64x2 8 -0.0079 0.0079
        const K2: F64x2 = F64x2::from_bits(0xBFE0000000000000, 0x3BEF620EAB9F01C8); // -4.999999999999999999468350684492e-1
        const K3: F64x2 = F64x2::from_bits(0x3FD5555555555555, 0x3C7520A75110AF2D); // 3.333333333333333331548473888204e-1
        const K4: F64x2 = F64x2::from_bits(0xBFD00000000000CD, 0x3C0D1BE0AFE27FF5); // -2.500000000000113795887528479459e-1
        const K5: F64x2 = F64x2::from_bits(0x3FC9999999999CC1, 0x3C6F16A00C4DE5CE); // 2.000000000000224233341875565267e-1
        const K6: F64x2 = F64x2::from_bits(0xBFC5555553EC31BC, 0x3C651D9EA2A5981B); // -1.666666660097585370166927936244e-1
        const K7: F64x2 = F64x2::from_bits(0x3FC24924904E0DE7, 0x3C61B0F5CC8E733D); // 1.428571419347541636820028854137e-1
        const K8: F64x2 = F64x2::from_bits(0xBFC00075C6134032, 0x3C5FCD3B2BFB0C68); // -1.250140397228292552397705588219e-1
        const K9: F64x2 = F64x2::from_bits(0x3FBC72C9FEF9C824, 0xBC523EFFDD8E4C78); // 1.111265418528835607196006034900e-1

        let (k, lo, ln_hi) = log_core_f64(x, edelta);
        let lo2 = lo.square();
        let ln_lo = lo + horner!(lo2, lo, [K2, K3, K4, K5, K6, K7, K8, K9]);

        // ln(x) = ln(2^k * m) = k * ln(2) + ln(m)
        k * LN_2 + (ln_hi + ln_lo)
    }

    #[inline]
    fn exp(x: F64x2) -> f64 {
        // x = k*ln(2) + r
        // k is an integer
        // |r| <= 0.5*ln(2)
        let (kf, k) = round_fi(x.hi() * LOG2_E);
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

        // exp(x) = exp(r) * 2^k = (t + 1) * 2^k
        (t + 1.0).scalbn_to_f64(k as i32)
    }

    // |z| = |x|^y = exp(y * ln(|x|))
    let lnx = ln(x.abs(), xedelta);
    let ylnx = y * lnx.hi();
    let absz = if ylnx.exponent() >= 10 {
        if ylnx.is_sign_negative() {
            0.0
        } else {
            f64::INFINITY
        }
    } else {
        exp(y * lnx)
    };

    absz.set_sign(sign)
}
