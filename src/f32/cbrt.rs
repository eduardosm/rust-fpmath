//! Cube root for `f32`.
//!
//! It is evaluated with plain `f64` arithmetic. The argument is split as
//! `|x| = 8^kdiv3 * 2^kmod3 * r0`, with an integer `kdiv3`, `0 <= kmod3 <= 2`
//! and `1 <= r0 < 2`, so
//!
//! `cbrt(|x|) = 2^kdiv3 * cbrt(2^kmod3) * cbrt(r0)`
//!
//! where `cbrt(2^kmod3)` is 1, `cbrt(2)` or `cbrt(4)`, and `cbrt(r0)` is
//! approximated with a polynomial, with a relative error of about 2^-10.5.
//! That approximation of `cbrt(|x|)` is refined with one Halley iteration,
//! which roughly cubes the relative error, so it is about 2^-32 before the
//! final rounding.

use crate::traits::Float as _;

// GENERATE: consts f64 CBRT_2 CBRT_4
const CBRT_2: f64 = f64::from_bits(0x3FF428A2F98D728B); // 1.2599210498948732e0
const CBRT_4: f64 = f64::from_bits(0x3FF965FEA53D6E3D); // 1.5874010519681996e0

// [1, cbrt(2), cbrt(4)]
static CBRT_SCALE: [f64; 3] = [1.0, CBRT_2, CBRT_4];

impl crate::generic::Cbrt for f32 {
    #[inline]
    fn cbrt_finite(x: Self) -> Self {
        let absx = f64::from(x.abs());

        // Split |x| = 2^k * r0 such as
        // * k is an integer
        // * 1 <= r0 < 2
        let k = absx.exponent();
        let r0 = absx.set_exp(0);

        // k = 3 * kdiv3 + kmod3, with 0 <= kmod3 <= 2
        let kp = (k + 153) as u16;
        let kmod3 = kp % 3;
        let kdiv3 = (kp / 3) as i16 - 51;

        // cbrt(|x|) = cbrt(r0 * 2^kmod3) * 2^kdiv3
        //           = cbrt(r0) * cbrt(2^kmod3) * 2^kdiv3
        // s = cbrt(2^kmod3) * 2^kdiv3
        let s = f64::from_bits(
            CBRT_SCALE[usize::from(kmod3)]
                .to_bits()
                .wrapping_add((kdiv3 as u64) << f64::MANT_BITS),
        );

        // y0 ~= cbrt(|x|)
        let y0 = {
            // GENERATE: cbrt_poly f64 3
            const K0: f64 = f64::from_bits(0x3FE3EE5A15141A55); // 6.22845688981054e-1
            const K1: f64 = f64::from_bits(0x3FDC0069D19700AE); // 4.375252291465682e-1
            const K2: f64 = f64::from_bits(0xBFAE8D1148F1473C); // -5.967000976001066e-2

            (K0 + horner!(r0, r0, [K1, K2])) * s
        };

        // Refine y0 with a Halley iteration:
        // y1 = y0 * (y0^3 + 2 * |x|) / (2 * y0^3 + |x|)
        let y03 = y0 * y0 * y0;
        let y1 = y0 * (y03 + 2.0 * absx) / (2.0 * y03 + absx);

        (y1 as f32).copysign(x)
    }
}
