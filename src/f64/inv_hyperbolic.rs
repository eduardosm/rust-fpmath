use super::f64x2::F64x2;
use super::log_core::log_core_f64x2;
use crate::traits::Float as _;

// GENERATE: consts F64x2 LN_2
const LN_2: F64x2 = F64x2::from_bits(0x3FE62E42FEFA39EF, 0x3C7ABC9E3B39803F); // 6.931471805599453094172321214582e-1

impl crate::generic::InvHyperbolic for f64 {
    fn asinh_finite(x: Self) -> Self {
        #[inline]
        fn ln(x: F64x2, edelta: i16) -> F64x2 {
            // GENERATE: ln_1p_poly F64x2 7 -0.0079 0.0079
            const K2: F64x2 = F64x2::from_bits(0xBFE0000000000000, 0x3BF48AE7307F28E9); // -4.999999999999999999303989972132e-1
            const K3: F64x2 = F64x2::from_bits(0x3FD55555555555AF, 0xBC785DAF7B868A9F); // 3.333333333333382896991266095150e-1
            const K4: F64x2 = F64x2::from_bits(0xBFD0000000000108, 0x3C5423562429FB02); // -2.500000000000146505771851914311e-1
            const K5: F64x2 = F64x2::from_bits(0x3FC999999890A76E, 0x3C6199AB1DC1AA7B); // 1.999999995180660815331568135373e-1
            const K6: F64x2 = F64x2::from_bits(0xBFC5555553A9951F, 0x3C647B49A8E55C3F); // -1.666666658885924323645098523820e-1
            const K7: F64x2 = F64x2::from_bits(0x3FC24994511B3AC0, 0xBC32B4480CCE6FB6); // 1.428704639460729442909190302719e-1
            const K8: F64x2 = F64x2::from_bits(0xBFC0007FC9817223, 0x3C670499BF4D9226); // -1.250152334131523113132157023371e-1

            let (k, lo, ln_hi) = log_core_f64x2(x, edelta);
            let lo2 = lo.square();
            let ln_lo = lo + horner!(lo2, lo, [K2, K3, K4, K5, K6, K7, K8]);

            // ln(x) = ln(2^k * m) = k * ln(2) + ln(m)
            k * LN_2 + (ln_hi + ln_lo)
        }

        #[inline]
        fn asinh_small(x: f64) -> F64x2 {
            // GENERATE: asinh_poly F64x2 4 -0.016 0.016
            const K3: F64x2 = F64x2::from_bits(0xBFC5555555555555, 0xBC478427F8538666); // -1.666666666666666599644434242178e-1
            const K5: F64x2 = F64x2::from_bits(0x3FB333333332EB99, 0x3C555B716C6CD57E); // 7.499999999974562200502960908291e-2
            const K7: F64x2 = F64x2::from_bits(0xBFA6DB6D9D13FA9E, 0xBC17A84064749E00); // -4.464285414177492116037335393196e-2
            const K9: F64x2 = F64x2::from_bits(0x3F9F18C94F9C2E6F, 0x3C38C488ABB26AFA); // 3.036799000169171157055376753311e-2

            let x2 = F64x2::square1(x);
            let x3 = x2 * x;
            x + horner!(x3, x2, [K3, K5, K7, K9])
        }

        if x.exponent() <= -7 {
            asinh_small(x).to_f64()
        } else if x.exponent() > 500 {
            // x is very large, avoid overflow when squaring x
            // |x| + sqrt(x^2 + 1) ~= 2 * |x|
            ln(F64x2::new1(x.abs()), 1).to_f64().copysign(x)
        } else {
            // t1 = |x| + sqrt(x^2 + 1)
            let t1 = x.abs() + (F64x2::square1(x) + 1.0).sqrt();

            // t2 = |asinh(x)| = ln(|x| + sqrt(x^2 + 1)) = ln(t1)
            let t2 = ln(t1, 0);

            // asinh(x) = |asinh(x)| * sgn(x)
            t2.to_f64().copysign(x)
        }
    }

    fn acosh_finite(x: Self) -> Self {
        #[inline]
        fn ln(x: F64x2, edelta: i16) -> F64x2 {
            // GENERATE: ln_1p_poly F64x2 7 -0.0079 0.0079
            const K2: F64x2 = F64x2::from_bits(0xBFE0000000000000, 0x3BF48AE7307F28E9); // -4.999999999999999999303989972132e-1
            const K3: F64x2 = F64x2::from_bits(0x3FD55555555555AF, 0xBC785DAF7B868A9F); // 3.333333333333382896991266095150e-1
            const K4: F64x2 = F64x2::from_bits(0xBFD0000000000108, 0x3C5423562429FB02); // -2.500000000000146505771851914311e-1
            const K5: F64x2 = F64x2::from_bits(0x3FC999999890A76E, 0x3C6199AB1DC1AA7B); // 1.999999995180660815331568135373e-1
            const K6: F64x2 = F64x2::from_bits(0xBFC5555553A9951F, 0x3C647B49A8E55C3F); // -1.666666658885924323645098523820e-1
            const K7: F64x2 = F64x2::from_bits(0x3FC24994511B3AC0, 0xBC32B4480CCE6FB6); // 1.428704639460729442909190302719e-1
            const K8: F64x2 = F64x2::from_bits(0xBFC0007FC9817223, 0x3C670499BF4D9226); // -1.250152334131523113132157023371e-1

            let (k, lo, ln_hi) = log_core_f64x2(x, edelta);
            let lo2 = lo.square();
            let ln_lo = lo + horner!(lo2, lo, [K2, K3, K4, K5, K6, K7, K8]);

            // ln(x) = ln(2^k * m) = k * ln(2) + ln(m)
            k * LN_2 + (ln_hi + ln_lo)
        }

        if x.exponent() > 500 {
            // x is very large, avoid overflow when squaring x
            // x + sqrt(x^2 - 1) ~= 2 * x
            ln(F64x2::new1(x), 1).to_f64()
        } else {
            // t1 = x + sqrt(x^2 - 1)
            let t1 = x + (F64x2::square1(x) - 1.0).sqrt();

            // acosh(x) = ln(x + sqrt(x^2 - 1)) = ln(t1)
            ln(t1, 0).to_f64()
        }
    }

    fn atanh_finite(x: Self) -> Self {
        #[inline]
        fn ln(x: F64x2) -> F64x2 {
            // GENERATE: ln_1p_poly F64x2 7 -0.0079 0.0079
            const K2: F64x2 = F64x2::from_bits(0xBFE0000000000000, 0x3BF48AE7307F28E9); // -4.999999999999999999303989972132e-1
            const K3: F64x2 = F64x2::from_bits(0x3FD55555555555AF, 0xBC785DAF7B868A9F); // 3.333333333333382896991266095150e-1
            const K4: F64x2 = F64x2::from_bits(0xBFD0000000000108, 0x3C5423562429FB02); // -2.500000000000146505771851914311e-1
            const K5: F64x2 = F64x2::from_bits(0x3FC999999890A76E, 0x3C6199AB1DC1AA7B); // 1.999999995180660815331568135373e-1
            const K6: F64x2 = F64x2::from_bits(0xBFC5555553A9951F, 0x3C647B49A8E55C3F); // -1.666666658885924323645098523820e-1
            const K7: F64x2 = F64x2::from_bits(0x3FC24994511B3AC0, 0xBC32B4480CCE6FB6); // 1.428704639460729442909190302719e-1
            const K8: F64x2 = F64x2::from_bits(0xBFC0007FC9817223, 0x3C670499BF4D9226); // -1.250152334131523113132157023371e-1

            let (k, lo, ln_hi) = log_core_f64x2(x, 0);
            let lo2 = lo.square();
            let ln_lo = lo + horner!(lo2, lo, [K2, K3, K4, K5, K6, K7, K8]);

            // ln(x) = ln(2^k * m) = k * ln(2) + ln(m)
            k * LN_2 + (ln_hi + ln_lo)
        }

        #[inline]
        fn ln_1p_small(x: F64x2) -> F64x2 {
            // GENERATE: ln_1p_poly F64x2 4 -0.0002442 0.0002442
            const K2: F64x2 = F64x2::from_bits(0xBFDFFFFFFFFFFFFD, 0x3C7586BFDDFAC003); // -4.999999999999998147954008087190e-1
            const K3: F64x2 = F64x2::from_bits(0x3FD555555555554F, 0x3C75FC1B3C3ACA71); // 3.333333333333330008314750331139e-1
            const K4: F64x2 = F64x2::from_bits(0xBFD000000D56ECCA, 0xBC5E3006F20FB4A1); // -2.500000124234246289223380012117e-1
            const K5: F64x2 = F64x2::from_bits(0x3FC99999B7DE5145, 0x3C1C03E76F6D5D76); // 2.000000140948349495650749156055e-1

            let x2 = x.square();
            x + horner!(x2, x, [K2, K3, K4, K5])
        }

        let absx = F64x2::new1(x.abs());

        if x.exponent() <= -14 {
            // t1 = (1 + |x|) / (1 - |x|) - 1 = 2 * |x| / (1 - |x|)
            let t1 = absx.twice() / (1.0 - absx);

            // atanh(x) = 0.5 * ln((1 + |x|) / (1 - |x|)) * sgn(x)
            //          = 0.5 * ln(t1 + 1) * sgn(x)
            ln_1p_small(t1).scalbn_to_f64(-1).copysign(x)
        } else {
            // t1 = (1 + |x|) / (1 - |x|)
            let t1 = (1.0 + absx) / (1.0 - absx);

            // atanh(x) = 0.5 * ln((1 + |x|) / (1 - |x|)) * sgn(x)
            //          = 0.5 * ln(t1) * sgn(x)
            ln(t1).halve().to_f64().copysign(x)
        }
    }
}
