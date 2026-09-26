use super::f64x2::F64x2;
use super::log_core::{log_core_f64, log_core_f64x2};
use crate::traits::Float as _;

// GENERATE: consts F64x2 LN_2 LOG2_E LOG10_E LOG10_2
const LN_2: F64x2 = F64x2::from_bits(0x3FE62E42FEFA39EF, 0x3C7ABC9E3B39803F); // 6.931471805599453094172321214582e-1
const LOG2_E: F64x2 = F64x2::from_bits(0x3FF71547652B82FE, 0x3C7777D0FFDA0D24); // 1.442695040888963407359924681002e0
const LOG10_E: F64x2 = F64x2::from_bits(0x3FDBCB7B1526E50E, 0x3C695355BAAAFAD3); // 4.342944819032518276511289189166e-1
const LOG10_2: F64x2 = F64x2::from_bits(0x3FD34413509F79FF, 0xBC49DC1DA994FD21); // 3.010299956639811952137388947245e-1

impl crate::generic::Log for f64 {
    fn ln_finite(x: Self, edelta: i16) -> Self {
        // GENERATE: ln_1p_poly F64x2 7 -0.0079 0.0079
        const K2: F64x2 = F64x2::from_bits(0xBFE0000000000000, 0x3BF48AE7307F28E9); // -4.999999999999999999303989972132e-1
        const K3: F64x2 = F64x2::from_bits(0x3FD55555555555AF, 0xBC785DAF7B868A9F); // 3.333333333333382896991266095150e-1
        const K4: F64x2 = F64x2::from_bits(0xBFD0000000000108, 0x3C5423562429FB02); // -2.500000000000146505771851914311e-1
        const K5: F64x2 = F64x2::from_bits(0x3FC999999890A76E, 0x3C6199AB1DC1AA7B); // 1.999999995180660815331568135373e-1
        const K6: F64x2 = F64x2::from_bits(0xBFC5555553A9951F, 0x3C647B49A8E55C3F); // -1.666666658885924323645098523820e-1
        const K7: F64x2 = F64x2::from_bits(0x3FC24994511B3AC0, 0xBC32B4480CCE6FB6); // 1.428704639460729442909190302719e-1
        const K8: F64x2 = F64x2::from_bits(0xBFC0007FC9817223, 0x3C670499BF4D9226); // -1.250152334131523113132157023371e-1

        // Split x * 2^edelta = 2^k * hi * (1 + lo)
        let (k, lo, ln_hi) = log_core_f64(x, edelta);

        // ln_lo = ln(1 + lo)
        let lo2 = lo.square();
        let ln_lo = lo + horner!(lo2, lo, [K2, K3, K4, K5, K6, K7, K8]);

        // x * 2^edelta = 2^k * hi * (1 + lo)
        // ln(x * 2^edelta) = k * ln(2) + ln(hi) + ln(1 + lo)
        let r = LN_2 * k + (ln_hi + ln_lo);
        r.to_f64()
    }

    fn ln_1p_finite(x: Self) -> Self {
        if x.exponent() <= -10 {
            // GENERATE: ln_1p_poly F64x2 5 -0.00196 0.00196
            const K2: F64x2 = F64x2::from_bits(0xBFE0000000000000, 0xBC32C58EA3DA65FF); // -5.000000000000000010176163169742e-1
            const K3: F64x2 = F64x2::from_bits(0x3FD555555555156D, 0xBC7CEE573A1E6D2F); // 3.333333333324251273018144086490e-1
            const K4: F64x2 = F64x2::from_bits(0xBFCFFFFFFFFEE9BC, 0xBC46E4E38CF29F6F); // -2.499999999980228062976280450238e-1
            const K5: F64x2 = F64x2::from_bits(0x3FC999A0027DBBEB, 0xBC660D15A9E6E191); // 2.000007640994846029717063025414e-1
            const K6: F64x2 = F64x2::from_bits(0xBFC5555D19BF2D98, 0x3C53B6A205DE8C33); // -1.666675925938363402309113121331e-1

            let x2 = F64x2::square1(x);
            let r = x + horner!(x2, x, [K2, K3, K4, K5, K6]);
            r.to_f64()
        } else {
            // GENERATE: ln_1p_poly F64x2 7 -0.0079 0.0079
            const K2: F64x2 = F64x2::from_bits(0xBFE0000000000000, 0x3BF48AE7307F28E9); // -4.999999999999999999303989972132e-1
            const K3: F64x2 = F64x2::from_bits(0x3FD55555555555AF, 0xBC785DAF7B868A9F); // 3.333333333333382896991266095150e-1
            const K4: F64x2 = F64x2::from_bits(0xBFD0000000000108, 0x3C5423562429FB02); // -2.500000000000146505771851914311e-1
            const K5: F64x2 = F64x2::from_bits(0x3FC999999890A76E, 0x3C6199AB1DC1AA7B); // 1.999999995180660815331568135373e-1
            const K6: F64x2 = F64x2::from_bits(0xBFC5555553A9951F, 0x3C647B49A8E55C3F); // -1.666666658885924323645098523820e-1
            const K7: F64x2 = F64x2::from_bits(0x3FC24994511B3AC0, 0xBC32B4480CCE6FB6); // 1.428704639460729442909190302719e-1
            const K8: F64x2 = F64x2::from_bits(0xBFC0007FC9817223, 0x3C670499BF4D9226); // -1.250152334131523113132157023371e-1

            let xp1 = F64x2::add11(x, 1.0);

            // Split x + 1 = 2^k * hi * (1 + lo)
            let (k, lo, ln_hi) = log_core_f64x2(xp1, 0);

            // ln_lo = ln(1 + lo)
            let lo2 = lo.square();
            let ln_lo = lo + horner!(lo2, lo, [K2, K3, K4, K5, K6, K7, K8]);

            // x * 2^edelta = 2^k * hi * (1 + lo)
            // ln(x * 2^edelta) = k * ln(2) + ln(hi) + ln(1 + lo)
            let r = LN_2 * k + (ln_hi + ln_lo);
            r.to_f64()
        }
    }

    fn log2_finite(x: Self, edelta: i16) -> Self {
        // GENERATE: ln_1p_poly F64x2 7 -0.0079 0.0079
        const K2: F64x2 = F64x2::from_bits(0xBFE0000000000000, 0x3BF48AE7307F28E9); // -4.999999999999999999303989972132e-1
        const K3: F64x2 = F64x2::from_bits(0x3FD55555555555AF, 0xBC785DAF7B868A9F); // 3.333333333333382896991266095150e-1
        const K4: F64x2 = F64x2::from_bits(0xBFD0000000000108, 0x3C5423562429FB02); // -2.500000000000146505771851914311e-1
        const K5: F64x2 = F64x2::from_bits(0x3FC999999890A76E, 0x3C6199AB1DC1AA7B); // 1.999999995180660815331568135373e-1
        const K6: F64x2 = F64x2::from_bits(0xBFC5555553A9951F, 0x3C647B49A8E55C3F); // -1.666666658885924323645098523820e-1
        const K7: F64x2 = F64x2::from_bits(0x3FC24994511B3AC0, 0xBC32B4480CCE6FB6); // 1.428704639460729442909190302719e-1
        const K8: F64x2 = F64x2::from_bits(0xBFC0007FC9817223, 0x3C670499BF4D9226); // -1.250152334131523113132157023371e-1

        // Split x * 2^edelta = 2^k * hi * (1 + lo)
        let (k, lo, ln_hi) = log_core_f64(x, edelta);

        // ln_lo = ln(1 + lo)
        let lo2 = lo.square();
        let ln_lo = lo + horner!(lo2, lo, [K2, K3, K4, K5, K6, K7, K8]);

        // x * 2^edelta = 2^k * hi * (1 + lo)
        // ln(x * 2^edelta) = k * ln(2) + ln(hi) + ln(1 + lo)
        let r = k + LOG2_E * (ln_hi + ln_lo);
        r.to_f64()
    }

    fn log10_finite(x: Self, edelta: i16) -> Self {
        // GENERATE: ln_1p_poly F64x2 7 -0.0079 0.0079
        const K2: F64x2 = F64x2::from_bits(0xBFE0000000000000, 0x3BF48AE7307F28E9); // -4.999999999999999999303989972132e-1
        const K3: F64x2 = F64x2::from_bits(0x3FD55555555555AF, 0xBC785DAF7B868A9F); // 3.333333333333382896991266095150e-1
        const K4: F64x2 = F64x2::from_bits(0xBFD0000000000108, 0x3C5423562429FB02); // -2.500000000000146505771851914311e-1
        const K5: F64x2 = F64x2::from_bits(0x3FC999999890A76E, 0x3C6199AB1DC1AA7B); // 1.999999995180660815331568135373e-1
        const K6: F64x2 = F64x2::from_bits(0xBFC5555553A9951F, 0x3C647B49A8E55C3F); // -1.666666658885924323645098523820e-1
        const K7: F64x2 = F64x2::from_bits(0x3FC24994511B3AC0, 0xBC32B4480CCE6FB6); // 1.428704639460729442909190302719e-1
        const K8: F64x2 = F64x2::from_bits(0xBFC0007FC9817223, 0x3C670499BF4D9226); // -1.250152334131523113132157023371e-1

        // Split x * 2^edelta = 2^k * hi * (1 + lo)
        let (k, lo, ln_hi) = log_core_f64(x, edelta);

        // ln_lo = ln(1 + lo)
        let lo2 = lo.square();
        let ln_lo = lo + horner!(lo2, lo, [K2, K3, K4, K5, K6, K7, K8]);

        // x * 2^edelta = 2^k * hi * (1 + lo)
        // ln(x * 2^edelta) = k * ln(2) + ln(hi) + ln(1 + lo)
        let r = LOG10_2 * k + LOG10_E * (ln_hi + ln_lo);
        r.to_f64()
    }
}
