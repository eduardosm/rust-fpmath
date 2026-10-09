//! Cube root for `f64`.
//!
//! The argument is split as `|x| = 8^kdiv3 * r`, with an integer `kdiv3` and
//! `1 <= r < 8`, so
//!
//! `cbrt(x) = sgn(x) * 2^kdiv3 * cbrt(r)`
//!
//! `cbrt(r)` is first approximated with a polynomial of the mantissa of `r`,
//! multiplied by 1, `cbrt(2)` or `cbrt(4)`, with a relative error of about
//! 2^-13.6. That approximation is rounded to a multiple of 2^-16 (`yc`), so
//! `yc^3` and `r - yc^3` are exact, and
//!
//! `cbrt(r) = yc * cbrt(1 + e)`, with `e = (r - yc^3) / yc^3`
//!
//! where `cbrt(1 + e) - 1` (about 2^-13.5 at most, in magnitude) is
//! approximated with another polynomial. The relative error before the final
//! rounding is about 2^-65.

use crate::traits::Float as _;

// GENERATE: consts f64 CBRT_2 CBRT_4
const CBRT_2: f64 = f64::from_bits(0x3FF428A2F98D728B); // 1.2599210498948732e0
const CBRT_4: f64 = f64::from_bits(0x3FF965FEA53D6E3D); // 1.5874010519681996e0

// [1, cbrt(2), cbrt(4)]
static CBRT_SCALE: [f64; 3] = [1.0, CBRT_2, CBRT_4];

impl crate::generic::Cbrt for f64 {
    #[inline]
    fn cbrt_finite(x: Self) -> Self {
        let (xn, edelta) = x.normalize_arg();

        // Split |x| = 2^k * r0 such as
        // * k is an integer
        // * 1 <= r0 < 2
        let k = xn.exponent() + edelta;
        let r0 = xn.abs().set_exp(0);

        // k = 3 * kdiv3 + kmod3, with 0 <= kmod3 <= 2
        let kp = (k + 1077) as u16;
        let kmod3 = kp % 3;
        let kdiv3 = (kp / 3) as i16 - 359;

        // cbrt(|x|) = cbrt(r) * 2^kdiv3, with 1 <= r = r0 * 2^kmod3 < 8
        let r = r0.set_exp(kmod3 as i16);

        // y ~= cbrt(r)
        let y = {
            // GENERATE: cbrt_poly f64 4
            const K0: f64 = f64::from_bits(0x3FE1B7EA3EA937F0); // 5.537005637977455e-1
            const K1: f64 = f64::from_bits(0x3FE2BA2945673001); // 5.852247577217896e-1
            const K2: f64 = f64::from_bits(0xBFC4B11B5561E2F4); // -1.6165486973820686e-1
            const K3: f64 = f64::from_bits(0x3F975B306FDD2341); // 2.2808796718240435e-2

            let t = K0 + r0 * (K1 + r0 * (K2 + r0 * K3));
            t * CBRT_SCALE[usize::from(kmod3)]
        };

        // Round y to a multiple of 2^-16.
        // Since 0.99 ~< y ~< 2.01, yc^2 and yc^3 are exact.
        let round_16 = f64::exp2i_fast(36);
        let yc = ((y + round_16).purify() - round_16).purify();

        // r = yc^3 * (1 + e)
        // cbrt(r) = yc * cbrt(1 + e)
        // `r - yc^3` is exact
        let yc3 = yc * yc * yc;
        let e = (r - yc3) / yc3;

        // c ~= cbrt(1 + e) - 1
        let c = {
            // GENERATE: cbrt_1p_poly f64 4 -0.000275 0.000275
            const K0: f64 = f64::from_bits(0x3FD5555555555555); // 3.333333333333333e-1
            const K1: f64 = f64::from_bits(0xBFBC71C71C71C718); // -1.1111111111111105e-1
            const K2: f64 = f64::from_bits(0x3FAF9ADD4FA7553F); // 6.172839734396484e-2
            const K3: f64 = f64::from_bits(0xBFA511E8E98E2CB0); // -4.1152266035248686e-2

            e * (K0 + e * (K1 + e * (K2 + e * K3)))
        };

        // cbrt(r) ~= yc * (1 + c)
        let t = yc + yc * c;

        // cbrt(x) = sgn(x) * cbrt(r) * 2^kdiv3
        (t * f64::exp2i_fast(kdiv3)).copysign(x)
    }
}
