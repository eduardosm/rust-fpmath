use super::log_core::{log_core_f32, log_core_f64};
use crate::traits::Float as _;

// GENERATE: consts f64 LN_2 LOG2_E LOG10_E LOG10_2
const LN_2: f64 = f64::from_bits(0x3FE62E42FEFA39EF); // 6.931471805599453e-1
const LOG2_E: f64 = f64::from_bits(0x3FF71547652B82FE); // 1.4426950408889634e0
const LOG10_E: f64 = f64::from_bits(0x3FDBCB7B1526E50E); // 4.342944819032518e-1
const LOG10_2: f64 = f64::from_bits(0x3FD34413509F79FF); // 3.010299956639812e-1

impl crate::generic::Log for f32 {
    fn ln_finite(x: Self, edelta: i16) -> Self {
        // GENERATE: ln_1p_poly f64 5 -0.032 0.032
        const K2: f64 = f64::from_bits(0xBFE000000002A753); // -5.000000000193076e-1
        const K3: f64 = f64::from_bits(0x3FD555550FF605CD); // 3.333332687253375e-1
        const K4: f64 = f64::from_bits(0xBFCFFFFED1D81A63); // -2.4999985929771915e-1
        const K5: f64 = f64::from_bits(0x3FC9A047AFADF4FF); // 2.0020385816670935e-1
        const K6: f64 = f64::from_bits(0xBFC55D6E20C696D6); // -1.6691376304986844e-1

        // Split x * 2^edelta = 2^k * hi * (1 + lo)
        let (k, lo, ln_hi) = log_core_f32(x, edelta);

        // ln_lo = ln(1 + lo)
        let lo2 = lo * lo;
        let ln_lo = lo + horner!(lo2, lo, [K2, K3, K4, K5, K6]);

        // x * 2^edelta = 2^k * hi * (1 + lo)
        // ln(x * 2^edelta) = k * ln(2) + ln(hi) + ln(1 + lo)
        let r = LN_2 * k + (ln_hi + ln_lo);
        r as f32
    }

    fn ln_1p_finite(x: Self) -> Self {
        if x.exponent() <= -10 {
            // GENERATE: ln_1p_poly f64 4 -0.00196 0.00196
            const K2: f64 = f64::from_bits(0xBFDFFFFFFFFFC9F1); // -4.999999999992318e-1
            const K3: f64 = f64::from_bits(0x3FD555555554F43F); // 3.3333333333195364e-1
            const K4: f64 = f64::from_bits(0xBFD000035B35DA72); // -2.5000080020200877e-1
            const K5: f64 = f64::from_bits(0x3FC999A1376977B1); // 2.0000090795195782e-1

            let x = f64::from(x);
            let x2 = x * x;
            let r = x + horner!(x2, x, [K2, K3, K4, K5]);
            r as f32
        } else {
            // GENERATE: ln_1p_poly f64 5 -0.032 0.032
            const K2: f64 = f64::from_bits(0xBFE000000002A753); // -5.000000000193076e-1
            const K3: f64 = f64::from_bits(0x3FD555550FF605CD); // 3.333332687253375e-1
            const K4: f64 = f64::from_bits(0xBFCFFFFED1D81A63); // -2.4999985929771915e-1
            const K5: f64 = f64::from_bits(0x3FC9A047AFADF4FF); // 2.0020385816670935e-1
            const K6: f64 = f64::from_bits(0xBFC55D6E20C696D6); // -1.6691376304986844e-1

            let xp1 = f64::from(x) + 1.0;

            // Split x + 1 = 2^k * hi * (1 + lo)
            let (k, lo, ln_hi) = log_core_f64(xp1, 0);

            // ln_lo = ln(1 + lo)
            let lo2 = lo * lo;
            let ln_lo = lo + horner!(lo2, lo, [K2, K3, K4, K5, K6]);

            // x * 2^edelta = 2^k * hi * (1 + lo)
            // ln(x * 2^edelta) = k * ln(2) + ln(hi) + ln(1 + lo)
            let r = LN_2 * k + (ln_hi + ln_lo);
            r as f32
        }
    }

    fn log2_finite(x: Self, edelta: i16) -> Self {
        // GENERATE: ln_1p_poly f64 5 -0.032 0.032
        const K2: f64 = f64::from_bits(0xBFE000000002A753); // -5.000000000193076e-1
        const K3: f64 = f64::from_bits(0x3FD555550FF605CD); // 3.333332687253375e-1
        const K4: f64 = f64::from_bits(0xBFCFFFFED1D81A63); // -2.4999985929771915e-1
        const K5: f64 = f64::from_bits(0x3FC9A047AFADF4FF); // 2.0020385816670935e-1
        const K6: f64 = f64::from_bits(0xBFC55D6E20C696D6); // -1.6691376304986844e-1

        // Split x * 2^edelta = 2^k * hi * (1 + lo)
        let (k, lo, ln_hi) = log_core_f32(x, edelta);

        // ln_lo = ln(1 + lo)
        let lo2 = lo * lo;
        let ln_lo = lo + horner!(lo2, lo, [K2, K3, K4, K5, K6]);

        // x * 2^edelta = 2^k * hi * (1 + lo)
        // ln(x * 2^edelta) = k * ln(2) + ln(hi) + ln(1 + lo)
        let r = k + LOG2_E * (ln_hi + ln_lo);
        r as f32
    }

    fn log10_finite(x: Self, edelta: i16) -> Self {
        // GENERATE: ln_1p_poly f64 5 -0.032 0.032
        const K2: f64 = f64::from_bits(0xBFE000000002A753); // -5.000000000193076e-1
        const K3: f64 = f64::from_bits(0x3FD555550FF605CD); // 3.333332687253375e-1
        const K4: f64 = f64::from_bits(0xBFCFFFFED1D81A63); // -2.4999985929771915e-1
        const K5: f64 = f64::from_bits(0x3FC9A047AFADF4FF); // 2.0020385816670935e-1
        const K6: f64 = f64::from_bits(0xBFC55D6E20C696D6); // -1.6691376304986844e-1

        // Split x * 2^edelta = 2^k * hi * (1 + lo)
        let (k, lo, ln_hi) = log_core_f32(x, edelta);

        // ln_lo = ln(1 + lo)
        let lo2 = lo * lo;
        let ln_lo = lo + horner!(lo2, lo, [K2, K3, K4, K5, K6]);

        // x * 2^edelta = 2^k * hi * (1 + lo)
        // ln(x * 2^edelta) = k * ln(2) + ln(hi) + ln(1 + lo)
        let r = LOG10_2 * k + LOG10_E * (ln_hi + ln_lo);
        r as f32
    }
}
