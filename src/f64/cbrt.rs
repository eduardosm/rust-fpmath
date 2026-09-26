use super::f64x2::F64x2;
use crate::traits::Float as _;

// GENERATE: consts f64 INV_CBRT_2 INV_CBRT_4
const INV_CBRT_2: f64 = f64::from_bits(0x3FE965FEA53D6E3D); // 7.937005259840998e-1
const INV_CBRT_4: f64 = f64::from_bits(0x3FE428A2F98D728B); // 6.299605249474366e-1

impl crate::generic::Cbrt for f64 {
    fn cbrt_finite(x: Self, edelta: Self::Exp) -> Self {
        // Split |x| * 2^edelta = 2^k * r such as
        // * k is an integer
        // * 1 <= r < 8
        let k = x.exponent() + edelta;
        // 0 <= kmod3 <= 2
        let kmod3 = (((k + 1077) as u16) % 3) as i16;

        // 1 <= r < 8
        let r = x.abs().set_exp(kmod3);

        // t1 ~= cbrt(1 / r) with a polynomial approximation
        // The polynomial approximation is for range [1, 2]
        // cbrt(1 / r) = cbrt(1 / r') * cbrt(1 / 2^kmod3)
        let t0 = {
            // GENERATE: inv_cbrt_poly f64 5
            const K0: f64 = f64::from_bits(0x3FFAC22A09449B4C); // 1.6724033700972596e0
            const K1: f64 = f64::from_bits(0xBFF2DFACC184CCFE); // -1.1796081122709547e0
            const K2: f64 = f64::from_bits(0x3FE67A424D6E9E1B); // 7.024241936059669e-1
            const K3: f64 = f64::from_bits(0xBFCCB7D044D303A5); // -2.2435954437790176e-1
            const K4: f64 = f64::from_bits(0x3F9DCE72DC97E5B4); // 2.9107851709317692e-2

            let r0 = r.set_exp(0);
            K0 + horner!(r0, r0, [K1, K2, K3, K4])
        };
        let t1 = match kmod3 {
            0 => t0,
            1 => t0 * INV_CBRT_2,
            2 => t0 * INV_CBRT_4,
            _ => unreachable!(),
        };

        // ti ~= cbrt(1 / r)
        // initially, ti = t1
        let ti = t1;

        // refine ti with Newton iterations
        // for each iteration: ti = ti - (1 / 3) * (r * ti^4 - ti)
        let frac_1_3 = F64x2::div11(1.0, 3.0);
        let ti = ti - frac_1_3 * (r * F64x2::square1(ti).square() - ti);
        let ti = ti - frac_1_3 * (r * ti.square().square() - ti);
        let ti = ti - frac_1_3 * (r * ti.square().square() - ti);

        // t2 = cbrt(r) = r * cbrt(1 / r)^2 = ti^2 * r
        let t2 = ti.square() * r;

        // cbrt(x) = sgn(x) * cbrt(r) * 2^((k - mod(k, 3)) / 3)
        //         = sgn(x) * t2 * 2^((k - mod(k, 3)) / 3)
        t2.scalbn_to_f64(((k - kmod3) / 3).into()).copysign(x)
    }
}
