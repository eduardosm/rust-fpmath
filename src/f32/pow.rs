use super::log_core::log_core_f32;
use crate::generic::{round_fi, scalbn};
use crate::traits::Float as _;

// GENERATE: consts f64 LN_2 LOG2_E
const LN_2: f64 = f64::from_bits(0x3FE62E42FEFA39EF); // 6.931471805599453e-1
const LOG2_E: f64 = f64::from_bits(0x3FF71547652B82FE); // 1.4426950408889634e0

impl crate::generic::Pow for f32 {
    fn pow_finite(x: Self, xedelta: Self::Exp, y: Self, sign: bool) -> Self {
        pow_core(x, xedelta, y.into(), sign)
    }

    fn powi_finite(x: Self, xedelta: Self::Exp, y: i32) -> Self {
        pow_core(x, xedelta, y.into(), x.is_sign_negative() && (y & 1) != 0)
    }
}

fn pow_core(x: f32, xedelta: i16, y: f64, sign: bool) -> f32 {
    #[inline]
    fn ln(x: f32, edelta: i16) -> f64 {
        // GENERATE: ln_1p_poly f64 6 -0.032 0.032
        const K2: f64 = f64::from_bits(0xBFE000000001FF45); // -5.000000000145312e-1
        const K3: f64 = f64::from_bits(0x3FD55555555FA9F6); // 3.3333333337091575e-1
        const K4: f64 = f64::from_bits(0xBFCFFFFF0A7CC9E0); // -2.4999988567431242e-1
        const K5: f64 = f64::from_bits(0x3FC999981E53B5B4); // 1.9999982338724254e-1
        const K6: f64 = f64::from_bits(0xBFC55CAB5870E137); // -1.6689054315953353e-1
        const K7: f64 = f64::from_bits(0x3FC25159F6260532); // 1.431076480767302e-1

        let (k, lo, ln_hi) = log_core_f32(x, edelta);
        let lo2 = lo * lo;
        let ln_lo = lo + horner!(lo2, lo, [K2, K3, K4, K5, K6, K7]);

        // ln(x) = ln(2^k * m) = k * ln(2) + ln(m)
        (k * LN_2 + ln_hi) + ln_lo
    }

    #[inline]
    fn exp(x: f64) -> f32 {
        if x.exponent() >= 7 {
            if x.is_sign_negative() {
                return 0.0;
            } else {
                return f32::INFINITY;
            }
        }

        // x = k*ln(2) + r
        // k is an integer
        // |r| <= 0.5*ln(2)
        let (kf, k) = round_fi(x * LOG2_E);
        let r = x - kf * LN_2;

        // exp(r) = exp(r / s)^s
        // s = 2^4

        let rs = r * (1.0 / 16.0);

        // GENERATE: exp_m1_poly f64 3 -0.022 0.022
        const K2: f64 = f64::from_bits(0x3FDFFFFFFFDAF480); // 4.9999999986523136e-1
        const K3: f64 = f64::from_bits(0x3FC555715CD78988); // 1.666700080232426e-1
        const K4: f64 = f64::from_bits(0x3FA555753B934B02); // 4.16676173423607e-2

        let rs2 = rs * rs;
        let mut t = rs + horner!(rs2, rs, [K2, K3, K4]);

        // exp(r) - 1 = (t + 1)^s - 1
        // (t + 1)^2 - 1 = t^2 + 2*t
        for _ in 0..4 {
            t = t * t + 2.0 * t;
        }

        // exp(x) = exp(r) * 2^k = (t + 1) * 2^k
        scalbn(t + 1.0, k as i32) as f32
    }

    // |z| = |x|^y = exp(y * ln(|x|))
    let absz = exp(y * ln(x.abs(), xedelta));
    absz.set_sign(sign)
}
