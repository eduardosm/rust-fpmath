use crate::generic::{round_fi, scalbn};
use crate::traits::Float as _;

// GENERATE: consts f32 LOG2_E LOG2_10
const LOG2_E: f32 = f32::from_bits(0x3FB8AA3B); // 1.442695e0
const LOG2_10: f32 = f32::from_bits(0x40549A78); // 3.321928e0

// GENERATE: consts f64 LN_2 LOG10_2 LN_10
const LN_2: f64 = f64::from_bits(0x3FE62E42FEFA39EF); // 6.931471805599453e-1
const LOG10_2: f64 = f64::from_bits(0x3FD34413509F79FF); // 3.010299956639812e-1
const LN_10: f64 = f64::from_bits(0x40026BB1BBB55516); // 2.302585092994046e0

impl crate::generic::Exp for f32 {
    fn exp_finite(x: Self) -> Self {
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
        let r = f64::from(x) - f64::from(kf) * LN_2;

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
        scalbn(t + 1.0, k) as f32
    }

    fn exp_m1_finite(x: Self) -> Self {
        if x.exponent() >= 7 {
            if x.is_sign_negative() {
                return -1.0;
            } else {
                return f32::INFINITY;
            }
        }

        // x = k*ln(2) + r
        // k is an integer
        // |r| <= 0.5*ln(2)
        let (kf, k) = round_fi(x * LOG2_E);
        let r = f64::from(x) - f64::from(kf) * LN_2;

        // exp(r) = exp(r / s)^s
        // s = 2^4

        let rs = r * (1.0 / 16.0);

        // GENERATE: exp_m1_poly f64 4 -0.022 0.022
        const K2: f64 = f64::from_bits(0x3FDFFFFFFFE414DC); // 4.99999999898433e-1
        const K3: f64 = f64::from_bits(0x3FC5555555449CEA); // 1.6666666663625246e-1
        const K4: f64 = f64::from_bits(0x3FA555718510AC6E); // 4.166750668840612e-2
        const K5: f64 = f64::from_bits(0x3F81112663DADA3F); // 8.333492204980362e-3

        let rs2 = rs * rs;
        let mut t = rs + horner!(rs2, rs, [K2, K3, K4, K5]);

        // exp(r) - 1 = (t + 1)^s - 1
        // (t + 1)^2 - 1 = t^2 + 2*t
        for _ in 0..4 {
            t = t * t + 2.0 * t;
        }

        // exp(x) - 1 = exp(r) * 2^k - 1 = t * 2^k + (2^k - 1)
        (scalbn(t, k) + (scalbn(1.0, k) - 1.0)) as f32
    }

    fn exp2_finite(x: Self) -> Self {
        if x.exponent() >= 8 {
            if x.is_sign_negative() {
                return 0.0;
            } else {
                return f32::INFINITY;
            }
        }

        // x = k + r*log2(e)
        // k is an integer
        // |r| <= 0.5*ln(2)
        let (kf, k) = round_fi(x);
        let r = (f64::from(x) - f64::from(kf)) * LN_2;

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
        scalbn(t + 1.0, k) as f32
    }

    fn exp10_finite(x: Self) -> Self {
        if x.exponent() >= 6 {
            if x.is_sign_negative() {
                return 0.0;
            } else {
                return f32::INFINITY;
            }
        }

        // x = k*log10(2) + r*log10(e)
        // k is an integer
        // |r| <= 0.5*ln(2)
        let (kf, k) = round_fi(x * LOG2_10);
        let r = (f64::from(x) - f64::from(kf) * LOG10_2) * LN_10;

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
        scalbn(t + 1.0, k) as f32
    }
}
