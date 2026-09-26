use super::log_core::log_core_f64;
use crate::traits::Float as _;

// GENERATE: consts f64 LN_2
const LN_2: f64 = f64::from_bits(0x3FE62E42FEFA39EF); // 6.931471805599453e-1

impl crate::generic::InvHyperbolic for f32 {
    fn asinh_finite(x: Self) -> Self {
        #[inline]
        fn ln(x: f64) -> f64 {
            // GENERATE: ln_1p_poly f64 5 -0.032 0.032
            const K2: f64 = f64::from_bits(0xBFE000000002A753); // -5.000000000193076e-1
            const K3: f64 = f64::from_bits(0x3FD555550FF605CD); // 3.333332687253375e-1
            const K4: f64 = f64::from_bits(0xBFCFFFFED1D81A63); // -2.4999985929771915e-1
            const K5: f64 = f64::from_bits(0x3FC9A047AFADF4FF); // 2.0020385816670935e-1
            const K6: f64 = f64::from_bits(0xBFC55D6E20C696D6); // -1.6691376304986844e-1

            let (k, lo, ln_hi) = log_core_f64(x, 0);
            let lo2 = lo * lo;
            let ln_lo = lo + horner!(lo2, lo, [K2, K3, K4, K5, K6]);

            // ln(x) = ln(2^k * m) = k * ln(2) + ln(m)
            k * LN_2 + (ln_hi + ln_lo)
        }

        if x.exponent() < -20 {
            // asinh(x) ~= x for small x
            return x;
        }

        let x = f64::from(x);

        // t1 = |x| + sqrt(x^2 + 1)
        let t1 = x.abs() + crate::f64::fast_sqrt(x * x + 1.0).0;

        // t2 = ln(t1) = |asinh(x)|

        // t2 = |asinh(x)| = ln(|x| + sqrt(x^2 + 1)) = ln(t1)
        let t2 = ln(t1);

        // asinh(x) = |asinh(x)| * sgn(x)
        t2.copysign(x) as f32
    }

    fn acosh_finite(x: Self) -> Self {
        #[inline]
        fn ln(x: f64) -> f64 {
            // GENERATE: ln_1p_poly f64 4 -0.032 0.032
            const K2: f64 = f64::from_bits(0xBFDFFFFFC5D4ED7F); // -4.999999458265946e-1
            const K3: f64 = f64::from_bits(0x3FD55554EC225CB3); // 3.3333323535903076e-1
            const K4: f64 = f64::from_bits(0xBFD0037D16C84A4D); // -2.50212929008886e-1
            const K5: f64 = f64::from_bits(0x3FC9A18851319890); // 2.0024208035028268e-1

            let (k, lo, ln_hi) = log_core_f64(x, 0);
            let lo2 = lo * lo;
            let ln_lo = lo + horner!(lo2, lo, [K2, K3, K4, K5]);

            // ln(x) = ln(2^k * m) = k * ln(2) + ln(m)
            k * LN_2 + (ln_hi + ln_lo)
        }

        let x = f64::from(x);

        // t1 = x + sqrt(x^2 - 1)
        let t1 = x + crate::f64::fast_sqrt(x * x - 1.0).0;

        // acosh(x) = ln(x + sqrt(x^2 - 1)) = ln(t1)
        ln(t1) as f32
    }

    fn atanh_finite(x: Self) -> Self {
        #[inline]
        fn ln(x: f64) -> f64 {
            // GENERATE: ln_1p_poly f64 5 -0.032 0.032
            const K2: f64 = f64::from_bits(0xBFE000000002A753); // -5.000000000193076e-1
            const K3: f64 = f64::from_bits(0x3FD555550FF605CD); // 3.333332687253375e-1
            const K4: f64 = f64::from_bits(0xBFCFFFFED1D81A63); // -2.4999985929771915e-1
            const K5: f64 = f64::from_bits(0x3FC9A047AFADF4FF); // 2.0020385816670935e-1
            const K6: f64 = f64::from_bits(0xBFC55D6E20C696D6); // -1.6691376304986844e-1

            let (k, lo, ln_hi) = log_core_f64(x, 0);
            let lo2 = lo * lo;
            let ln_lo = lo + horner!(lo2, lo, [K2, K3, K4, K5, K6]);

            // ln(x) = ln(2^k * m) = k * ln(2) + ln(m)
            k * LN_2 + (ln_hi + ln_lo)
        }

        if x.exponent() < -20 {
            // atanh(x) ~= x for small x
            return x;
        }

        let x = f64::from(x);
        let absx = x.abs();

        // t1 = (1 + |x|) / (1 - |x|)
        let t1 = (1.0 + absx) / (1.0 - absx);

        // atanh(x) =
        //          = 0.5 * ln((1 + |x|) / (1 - |x|)) * sgn(x)
        //          = 0.5 * ln(t1) * sgn(x)
        (0.5 * ln(t1)).copysign(x) as f32
    }
}
