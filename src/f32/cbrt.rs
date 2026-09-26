use crate::generic::scalbn_medium;
use crate::traits::Float as _;

// GENERATE: consts f32 INV_CBRT_2 INV_CBRT_4
const INV_CBRT_2: f32 = f32::from_bits(0x3F4B2FF5); // 7.937005e-1
const INV_CBRT_4: f32 = f32::from_bits(0x3F214518); // 6.2996054e-1

impl crate::generic::Cbrt for f32 {
    fn cbrt_finite(x: Self, edelta: Self::Exp) -> Self {
        // Split |x| * 2^edelta = 2^k * r such as
        // * k is an integer
        // * 1 <= r < 8
        let k = x.exponent() + edelta;
        // 0 <= kmod3 <= 2
        let kmod3 = (((k + 153) as u16) % 3) as i16;

        // 1 <= r < 8
        let r = x.abs().set_exp(kmod3);

        // t1 ~= cbrt(1 / r) with a polynomial approximation
        // The polynomial approximation is for range [1, 2]
        // cbrt(1 / r) = cbrt(1 / r') * cbrt(1 / 2^kmod3)
        let t0 = {
            // GENERATE: inv_cbrt_poly f32 3
            const K0: f32 = f32::from_bits(0x3FB21939); // 1.3913947e0
            const K1: f32 = f32::from_bits(0xBEF9C752); // -4.8784882e-1
            const K2: f32 = f32::from_bits(0x3DC257A9); // 9.489376e-2

            let r0 = r.set_exp(0);
            K0 + horner!(r0, r0, [K1, K2])
        };
        let t1 = match kmod3 {
            0 => t0,
            1 => t0 * INV_CBRT_2,
            2 => t0 * INV_CBRT_4,
            _ => unreachable!(),
        };

        let r = f64::from(r);

        // ti ~= cbrt(1 / r)
        // initially, ti = t1
        let mut ti = f64::from(t1);

        // refine ti with Newton iterations
        // for each iteration: ti = ti - (1 / 3) * (r * ti^4 - ti)
        let frac_1_3 = 1.0 / 3.0;
        ti = ti - frac_1_3 * (r * pow4(ti) - ti);
        ti = ti - frac_1_3 * (r * pow4(ti) - ti);
        ti = ti - frac_1_3 * (r * pow4(ti) - ti);

        // t2 = cbrt(r) = r * cbrt(1 / r)^2 = ti^2 * r
        let t2 = ti * ti * r;

        // y = cbrt(x)
        //   = sgn(x) * cbrt(r) * 2^((k - mod(k, 3)) / 3)
        //   = sgn(x) * t2 * 2^((k - mod(k, 3)) / 3)
        let y = scalbn_medium(t2, ((k - kmod3) / 3).into()).set_sign(x.is_sign_negative());

        y as f32
    }
}

#[inline]
fn pow4(x: f64) -> f64 {
    let x2 = x * x;
    x2 * x2
}
