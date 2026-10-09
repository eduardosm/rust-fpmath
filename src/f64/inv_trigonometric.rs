//! Inverse trigonometric functions for `f64`.
//!
//! # `atan` and `atan2`
//!
//! They are calculated from the arc tangent of a ratio `n / d`, with
//! `0 <= n <= d` (see `atan_parts`), plus a multiple of π/2:
//! * `atan(x)`: `|x|`, or `1 / |x|` when `|x| > 1`, with
//!   `atan(|x|) = π/2 - atan(1 / |x|)`.
//! * `atan2(y, x)`: `min(|x|, |y|) / max(|x|, |y|)`.
//!
//! The arc tangent is reduced as
//!
//! `atan(n / d) = atan(c) + atan(t)`, with `t = (n - c * d) / (d + c * n)`
//!
//! where `c = i / 64` is the closest to `n / d`, so `|t| <= 1/128`.
//! `atan(c)` is taken from `ATAN_TBL` as a double-word, `t` is calculated as
//! a double-word, and `atan(t) - t` is approximated with a short polynomial.
//! The leading parts are added exactly, so the relative error before the
//! final rounding is about 2^-65.
//!
//! # `asin` and `acos`
//!
//! For `|x| <= 1/2`, `asin(x) = x + x * S(x^2)`, where `S` is approximated
//! with a polynomial (see `asin_poly`), and `acos(x) = π/2 - asin(x)`.
//!
//! For `|x| > 1/2`, `asin(|x|) = π/2 - 2 * asin(y)`, `acos(|x|) = 2 * asin(y)`
//! and `acos(-|x|) = π - 2 * asin(y)`, with `y = sqrt((1 - |x|) / 2)`, which
//! is calculated as a double-word, and `asin(y) = y + y * S(y^2)`, where
//! `y^2 = (1 - |x|) / 2` is exact.
//!
//! The relative error before the final rounding is about 2^-63.
//!
//! # Degrees and half-turns
//!
//! The functions in degrees and half-turns multiply the result (as a
//! double-word) by `180 / π` or `1 / π`.

use super::f64x2::{F64x2, Split};
use super::{round_i32, sqrt_parts, square_parts};
use crate::traits::Float as _;

// GENERATE: consts F64x2 FRAC_PI_2 PI FRAC_180_PI FRAC_1_PI
const FRAC_PI_2: F64x2 = F64x2::from_bits(0x3FF921FB54442D18, 0x3C91A62633145C07); // 1.570796326794896619231321691640e0
const PI: F64x2 = F64x2::from_bits(0x400921FB54442D18, 0x3CA1A62633145C07); // 3.141592653589793238462643383280e0
const FRAC_180_PI: F64x2 = F64x2::from_bits(0x404CA5DC1A63C1F8, 0xBCE1E7AB456405F9); // 5.729577951308232087679815481411e1
const FRAC_1_PI: F64x2 = F64x2::from_bits(0x3FD45F306DC9C883, 0xBC76B01EC5417056); // 3.183098861837906715377675267450e-1

// GENERATE: atan_table F64x2 6
// ATAN_TBL[i] = atan(i / 64)
static ATAN_TBL: [F64x2; 65] = [
    F64x2::from_bits(0x0000000000000000, 0x0000000000000000), // 0
    F64x2::from_bits(0x3F8FFF555BBB729B, 0xBC2220C39D4DFF50), // 1.562372862047683080280152125657e-2
    F64x2::from_bits(0x3F9FFD55BBA97625, 0xBC35EC431444912C), // 3.123983343026827625371174489249e-2
    F64x2::from_bits(0x3FA7FB818430DA2A, 0xBC086EF8F794F105), // 4.684071291596965375222376000186e-2
    F64x2::from_bits(0x3FAFF55BB72CFDEA, 0xBC3C934D86D23F1D), // 6.241880999595734847397911298551e-2
    F64x2::from_bits(0x3FB3F59F0E7C559D, 0x3C5AC4CE285DF847), // 7.796663383154230656332864878103e-2
    F64x2::from_bits(0x3FB7EE182602F10F, 0xBC5CFB654C0C3D98), // 9.347678115858946350452719331206e-2
    F64x2::from_bits(0x3FBBE39EBE6F07C3, 0x3C5F7B8F29A05987), // 1.089419569898657998418608611626e-1
    F64x2::from_bits(0x3FBFD5BA9AAC2F6E, 0xBC4CD37686760C17), // 1.243549945467614350313548491639e-1
    F64x2::from_bits(0x3FC1E1FAFB043727, 0xBC4B485914DACF8C), // 1.397088742891636451833677767391e-1
    F64x2::from_bits(0x3FC3D6EEE8C6626C, 0x3C661A3B0CE9281B), // 1.549967419239409823037143749335e-1
    F64x2::from_bits(0x3FC5C9811E3EC26A, 0xBC5054AB2C010F3D), // 1.702119252854744044904966070998e-1
    F64x2::from_bits(0x3FC7B97B4BCE5B02, 0x3C5347B0B4F881CA), // 1.853479499956947648860259612285e-1
    F64x2::from_bits(0x3FC9A6A8E96C8626, 0x3C4CF601E7B4348E), // 2.003985538258785146539457850344e-1
    F64x2::from_bits(0x3FCB90D7529260A2, 0x3C217B10D2E0E5AB), // 2.153576996977380480244596271665e-1
    F64x2::from_bits(0x3FCD77D5DF205736, 0x3C6C648D1534597E), // 2.302195872768437302401709596798e-1
    F64x2::from_bits(0x3FCF5B75F92C80DD, 0x3C68AB6E3CF7AFBD), // 2.449786631268641541720824812113e-1
    F64x2::from_bits(0x3FD09DC597D86362, 0x3C762E47390CB865), // 2.596296294082575310299464431840e-1
    F64x2::from_bits(0x3FD18BF5A30BF178, 0x3C630CA4748B1BF9), // 2.741674511196587975993718983422e-1
    F64x2::from_bits(0x3FD278372057EF46, 0xBC7077CDD36DFC81), // 2.885873618940773956236114199582e-1
    F64x2::from_bits(0x3FD362773707EBCC, 0xBC6963A544B672D8), // 3.028848683749714055605560945056e-1
    F64x2::from_bits(0x3FD44AA436C2AF0A, 0xBC75D5E43C55B3BA), // 3.170557532091470098090155766745e-1
    F64x2::from_bits(0x3FD530AD9951CD4A, 0xBC62566480884082), // 3.310960767041320949443387877569e-1
    F64x2::from_bits(0x3FD614840309CFE2, 0xBC7A725715711F00), // 3.450021772071051088676812869001e-1
    F64x2::from_bits(0x3FD6F61941E4DEF1, 0xBC7C63AAE6F6E918), // 3.587706702705722203959200639265e-1
    F64x2::from_bits(0x3FD7D5604B63B3F7, 0x3C769C885C2B249A), // 3.723984466767542219236550382837e-1
    F64x2::from_bits(0x3FD8B24D394A1B25, 0x3C7B6D0BA3748FA8), // 3.858826693980737758976954846072e-1
    F64x2::from_bits(0x3FD98CD5454D6B18, 0x3C79E6C988FD0A77), // 3.992207695752525656147166961589e-1
    F64x2::from_bits(0x3FDA64EEC3CC23FD, 0xBC724DEC1B50B7FF), // 4.124104415973873068997912896671e-1
    F64x2::from_bits(0x3FDB3A911DA65C6C, 0x3C7AE187B1CA5040), // 4.254496373700422895422636051808e-1
    F64x2::from_bits(0x3FDC0DB4C94EC9F0, 0xBC7CC1CE70934C34), // 4.383365598579578054456160492148e-1
    F64x2::from_bits(0x3FDCDE53432C1351, 0xBC7A2CFA4418F1AD), // 4.510696559885234763756392572822e-1
    F64x2::from_bits(0x3FDDAC670561BB4F, 0x3C7A2B7F222F65E2), // 4.636476090008061162142562314612e-1
    F64x2::from_bits(0x3FDE77EB7F175A34, 0x3C70E53DC1BF3435), // 4.760693303227612340751004202615e-1
    F64x2::from_bits(0x3FDF40DD0B541418, 0xBC6A3992DC382A23), // 4.883339510564055238671649607471e-1
    F64x2::from_bits(0x3FE0039C73C1A40C, 0xBC8B32C949C9D593), // 5.004408131472941140300005149792e-1
    F64x2::from_bits(0x3FE0657E94DB30D0, 0xBC7D5B495F6349E6), // 5.123894603107377066666010205843e-1
    F64x2::from_bits(0x3FE0C6145B5B43DA, 0x3C5974FA13B5404F), // 5.241796287829132483216496175045e-1
    F64x2::from_bits(0x3FE1255D9BFBD2A9, 0xBC52BDAEE1C0EE35), // 5.358112379604637002690850687077e-1
    F64x2::from_bits(0x3FE1835A88BE7C13, 0x3C8C621CEC00C301), // 5.472843809874369739852207703128e-1
    F64x2::from_bits(0x3FE1E00BABDEFEB4, 0xBC5928DF287A668F), // 5.585993153435624359715082164017e-1
    F64x2::from_bits(0x3FE23B71E2CC9E6A, 0x3C6C421C9F38224E), // 5.697564534829784433238348916656e-1
    F64x2::from_bits(0x3FE2958E59308E31, 0xBC709E73B0C6C087), // 5.807563535676703992032744750015e-1
    F64x2::from_bits(0x3FE2EE628406CBCA, 0x3C8C5D5E9FF0CF8D), // 5.915997103351114331458526589590e-1
    F64x2::from_bits(0x3FE345F01CCE37BB, 0x3C81021137C71102), // 6.022873461349641816821226942042e-1
    F64x2::from_bits(0x3FE39C391CD4171A, 0xBC82304331D8BF46), // 6.128202021652413251433846354957e-1
    F64x2::from_bits(0x3FE3F13FB89E96F4, 0x3C7ECF8B492644F0), // 6.231993299340659309924753490604e-1
    F64x2::from_bits(0x3FE445065B795B56, 0xBC7F76D0163F79C8), // 6.334258829691445662686954830593e-1
    F64x2::from_bits(0x3FE4978FA3269EE1, 0x3C72419A87F2A458), // 6.435011087932843868028092287173e-1
    F64x2::from_bits(0x3FE4E8DE5BB6EC04, 0x3C84A33DBEB3796C), // 6.534263411807619628638934113116e-1
    F64x2::from_bits(0x3FE538F57B89061F, 0xBC81BB74ABDA520C), // 6.632029927060932553632543102383e-1
    F64x2::from_bits(0x3FE587D81F732FBB, 0xBC75E5C9D8C5A950), // 6.728325475937631893114013292615e-1
    F64x2::from_bits(0x3FE5D58987169B18, 0x3C60028E4BC5E7CA), // 6.823165548747480782564299817112e-1
    F64x2::from_bits(0x3FE6220D115D7B8E, 0xBC62B785350EE8C1), // 6.916566218531998629800663181104e-1
    F64x2::from_bits(0x3FE66D663923E087, 0xBC76EA6FEBE8BBBA), // 7.008544078844501724579512817868e-1
    F64x2::from_bits(0x3FE6B798920B3D99, 0xBC8A80386188C50E), // 7.099116184635248611916111509362e-1
    F64x2::from_bits(0x3FE700A7C5784634, 0xBC78C34D25AADEF6), // 7.188299996216245054170141515259e-1
    F64x2::from_bits(0x3FE748978FBA8E0F, 0x3C47B2A6165884A1), // 7.276113326265106787829526909500e-1
    F64x2::from_bits(0x3FE78F6BBD5D315E, 0x3C8406A089803740), // 7.362574289814281317428352710891e-1
    F64x2::from_bits(0x3FE7D528289FA093, 0x3C8560821E2F3AA9), // 7.447701257160751857639310909741e-1
    F64x2::from_bits(0x3FE819D0B7158A4D, 0xBC7BF76229D3B917), // 7.531512809621943895247393702690e-1
    F64x2::from_bits(0x3FE85D69576CC2C5, 0x3C66B66E7FC8B8C3), // 7.614027698055784264231855420876e-1
    F64x2::from_bits(0x3FE89FF5FF57F1F8, 0xBC855B9A5E177A1B), // 7.695264804056582604068200359857e-1
    F64x2::from_bits(0x3FE8E17AA99CC05E, 0xBC7EC182AB042F61), // 7.775243103733477667249308161224e-1
    F64x2::from_bits(0x3FE921FB54442D18, 0x3C81A62633145C07), // 7.853981633974483096156608458199e-1
];

impl crate::generic::InvTrigonometric for f64 {
    // GENERATE: consts f64 PI FRAC_PI_2
    const PI: f64 = f64::from_bits(0x400921FB54442D18); // 3.141592653589793e0
    const FRAC_PI_2: f64 = f64::from_bits(0x3FF921FB54442D18); // 1.5707963267948966e0

    #[inline]
    fn asin_finite(x: Self) -> Self {
        let (hi, lo) = asin_core(x);
        hi + lo
    }

    #[inline]
    fn acos_finite(x: Self) -> Self {
        let (hi, lo) = acos_core(x);
        hi + lo
    }

    #[inline]
    fn atan_finite(x: Self) -> Self {
        let (hi, lo) = atan_core(x);
        hi + lo
    }

    #[inline]
    fn atan2_finite(y: Self, x: Self) -> Self {
        let (hi, lo, edelta) = atan2_core(y, x);
        if edelta == 0 {
            hi + lo
        } else {
            F64x2::fast_add11(hi, lo).scalbn_to_f64(edelta)
        }
    }

    #[inline]
    fn asind_finite(x: Self) -> Self {
        if x.exponent() < -960 {
            mul_tiny(x, FRAC_180_PI)
        } else {
            let (hi, lo) = asin_core(x);
            mul_to_f64(hi, lo, FRAC_180_PI)
        }
    }

    #[inline]
    fn acosd_finite(x: Self) -> Self {
        let (hi, lo) = acos_core(x);
        mul_to_f64(hi, lo, FRAC_180_PI)
    }

    #[inline]
    fn atand_finite(x: Self) -> Self {
        if x.exponent() < -960 {
            mul_tiny(x, FRAC_180_PI)
        } else {
            let (hi, lo) = atan_core(x);
            mul_to_f64(hi, lo, FRAC_180_PI)
        }
    }

    #[inline]
    fn atan2d_finite(y: Self, x: Self) -> Self {
        let (hi, lo, edelta) = atan2_core(y, x);
        if edelta == 0 {
            mul_to_f64(hi, lo, FRAC_180_PI)
        } else {
            // `F64x2` can lose the sign of zero, `copysign` recovers it.
            (F64x2::fast_add11(hi, lo) * FRAC_180_PI)
                .scalbn_to_f64(edelta)
                .copysign(y)
        }
    }

    #[inline]
    fn asinpi_finite(x: Self) -> Self {
        if x.exponent() < -960 {
            mul_tiny(x, FRAC_1_PI)
        } else {
            let (hi, lo) = asin_core(x);
            mul_to_f64(hi, lo, FRAC_1_PI)
        }
    }

    #[inline]
    fn acospi_finite(x: Self) -> Self {
        let (hi, lo) = acos_core(x);
        mul_to_f64(hi, lo, FRAC_1_PI)
    }

    #[inline]
    fn atanpi_finite(x: Self) -> Self {
        if x.exponent() < -960 {
            mul_tiny(x, FRAC_1_PI)
        } else {
            let (hi, lo) = atan_core(x);
            mul_to_f64(hi, lo, FRAC_1_PI)
        }
    }

    #[inline]
    fn atan2pi_finite(y: Self, x: Self) -> Self {
        let (hi, lo, edelta) = atan2_core(y, x);
        if edelta == 0 {
            mul_to_f64(hi, lo, FRAC_1_PI)
        } else {
            // `F64x2` can lose the sign of zero, `copysign` recovers it.
            (F64x2::fast_add11(hi, lo) * FRAC_1_PI)
                .scalbn_to_f64(edelta)
                .copysign(y)
        }
    }
}

/// Returns `(hi + lo) * k` rounded to `f64`, where `|lo| <= 2^-14 * |hi|`.
#[inline]
fn mul_to_f64(hi: f64, lo: f64, k: F64x2) -> f64 {
    // (hi + lo) * k ~= (hi + lo) * k.hi() + hi * k.lo(). The error is about
    // 2^-65, mostly from the terms with `lo` (up to 2^-14 * |hi|), which are
    // rounded or omitted (`lo * k.lo()`). `k` is a constant at the call sites,
    // so the optimizer can fold `Split::new(k.hi())`.
    let (p_hi, p_lo) = Split::new(k.hi()).mul(hi, lo);
    p_hi + (p_lo + hi * k.lo())
}

/// Returns `x * k` rounded to `f64`, for `|x| < 2^-960`, where
/// `asin(x) ~= atan(x) ~= x`.
#[inline]
fn mul_tiny(x: f64, k: F64x2) -> f64 {
    // The product is scaled to avoid subnormals in `F64x2`, which break its
    // accuracy, and `copysign` recovers the sign of zero, which `F64x2` can
    // lose.
    (F64x2::new1(x * f64::exp2i_fast(106)) * k)
        .scalbn_to_f64(-106)
        .copysign(x)
}

/// Returns `(hi, lo)` such that `hi + lo ~= asin(x)`, for `|x| < 1`.
///
/// `|lo| <= 2^-14 * |hi|`.
#[inline]
fn asin_core(x: f64) -> (f64, f64) {
    let ax = x.abs();
    let (hi, lo) = if ax <= 0.5 {
        if ax < f64::exp2i_fast(-32) {
            // asin(x) = x + x^3 / 6 + ..., where the second term is less than
            // 2^-66.6 relative to the first one
            (ax, 0.0)
        } else {
            // asin(x) = x + p, where `|p| < 0.048 * |x|`
            let (p_hi, p_lo) = asin_poly_small(ax);
            let s = F64x2::fast_add11(ax, p_hi);
            (s.hi(), s.lo() + p_lo)
        }
    } else {
        // asin(|x|) = π/2 - 2 * asin(y) = π/2 - 2 * y - 2 * p, where
        // π/2 > 2 * y > 2 * p (y <= 0.5)
        let (y_hi, y_lo, p_hi, p_lo) = asin_half(ax);
        let s1 = F64x2::fast_add11(FRAC_PI_2.hi(), -2.0 * y_hi);
        let s2 = F64x2::fast_add11(s1.hi(), -2.0 * p_hi);
        let lo = (s1.lo() + s2.lo()) + (FRAC_PI_2.lo() - 2.0 * (y_lo + p_lo));
        (s2.hi(), lo)
    };
    if x.is_sign_negative() {
        (-hi, -lo)
    } else {
        (hi, lo)
    }
}

/// Returns `(hi, lo)` such that `hi + lo ~= acos(x)`, for `|x| < 1`.
///
/// `|lo| <= 2^-14 * |hi|`.
#[inline]
fn acos_core(x: f64) -> (f64, f64) {
    let ax = x.abs();
    if ax <= 0.5 {
        if ax < f64::exp2i_fast(-32) {
            // acos(x) = π/2 - x - x^3 / 6 - ..., where the third term is
            // negligible
            return (FRAC_PI_2.hi(), FRAC_PI_2.lo() - x);
        }
        // acos(x) = π/2 - asin(x) = π/2 - x - p, where π/2 > |x| and
        // |π/2 - x| > 1 > |p|
        let (p_hi, p_lo) = asin_poly_small(ax);
        let (p_hi, p_lo) = if x.is_sign_negative() {
            (-p_hi, -p_lo)
        } else {
            (p_hi, p_lo)
        };
        let s1 = F64x2::fast_add11(FRAC_PI_2.hi(), -x);
        let s2 = F64x2::fast_add11(s1.hi(), -p_hi);
        (s2.hi(), (s1.lo() + s2.lo()) + (FRAC_PI_2.lo() - p_lo))
    } else {
        let (y_hi, y_lo, p_hi, p_lo) = asin_half(ax);
        if x.is_sign_negative() {
            // acos(x) = π - 2 * asin(y) = π - 2 * y - 2 * p, where
            // π > 2 * y > 2 * p (y <= 0.5)
            let s1 = F64x2::fast_add11(PI.hi(), -2.0 * y_hi);
            let s2 = F64x2::fast_add11(s1.hi(), -2.0 * p_hi);
            let lo = (s1.lo() + s2.lo()) + (PI.lo() - 2.0 * (y_lo + p_lo));
            (s2.hi(), lo)
        } else {
            // acos(x) = 2 * asin(y) = 2 * y + 2 * p
            let s = F64x2::fast_add11(2.0 * y_hi, 2.0 * p_hi);
            (s.hi(), s.lo() + 2.0 * (y_lo + p_lo))
        }
    }
}

/// Returns `(p_hi, p_lo)` such that `p_hi + p_lo ~= asin(x) - x = x * S(x^2)`,
/// for `2^-32 <= x <= 1/2`, with an error of about 2^-64.5 relative to `x`.
#[inline]
fn asin_poly_small(x: f64) -> (f64, f64) {
    // u = x^2 as a normalized double-word
    let (u_hi, u_lo) = square_parts(x);
    let (s_hi, s_lo) = asin_poly(u_hi, u_lo);
    Split::new(x).mul(s_hi, s_lo)
}

/// Returns `(y_hi, y_lo, p_hi, p_lo)`, where `y_hi + y_lo ~= y` and
/// `p_hi + p_lo ~= asin(y) - y = y * S(y^2)`, with `y = sqrt((1 - x) / 2)`, for
/// `1/2 < x < 1`, with errors of about 2^-100 and 2^-64.5 relative to `y`.
#[inline]
fn asin_half(x: f64) -> (f64, f64, f64, f64) {
    // y^2 = z = (1 - x) / 2, exact (Sterbenz lemma), and z >= 2^-54 is normal
    let z = (1.0 - x) * 0.5;

    // y = sqrt(z) as a double-word
    let (y, y_lo) = sqrt_parts(z, 0.0);
    let ys = Split::new(y);

    // asin(y) - y = y * S(z), where `y_lo * S(z)` (less than 2^-54 * y) is
    // added to the low part
    let (s_hi, s_lo) = asin_poly(z, 0.0);
    let (p_hi, p_lo) = ys.mul(s_hi, s_lo);
    (y, y_lo, p_hi, p_lo + y_lo * s_hi)
}

/// Returns `(hi, lo)` such that `hi + lo ~= S(u) = asin(x) / x - 1`, where
/// `u = x^2 = u_hi + u_lo` (with `|u_lo| <= 2^-53 * u_hi`), for
/// `0 <= u <= 1/4`, with an absolute error of about 2^-64.5.
///
/// `S(u) < 0.048` and `|lo| <= 2^-23 * |hi|`.
#[inline]
fn asin_poly(u_hi: f64, u_lo: f64) -> (f64, f64) {
    // S(u) ~= u * (K1 + K2 * u + ... + K15 * u^14), with a relative error of
    // 2^-67.7 in `asin(x)`
    // GENERATE: asin_poly F64x2:3,f64 15 0.5
    const K1: F64x2 = F64x2::from_bits(0x3FC5555555555556, 0xBC6E79DD537C2A84); // 1.666666666666666719535110873307e-1
    const K2: F64x2 = F64x2::from_bits(0x3FB3333333333298, 0x3C45B4816E305F15); // 7.499999999999784852060387934620e-2
    const K3: F64x2 = F64x2::from_bits(0x3FA6DB6DB6DC1A96, 0xBC28E6FF551FD7C5); // 4.464285714316422392315172545692e-2
    const K4: f64 = f64::from_bits(0x3F9F1C71C6BAD0FA); // 3.0381944422246378e-2
    const K5: f64 = f64::from_bits(0x3F96E8BA3EF02307); // 2.2372160045082825e-2
    const K6: f64 = f64::from_bits(0x3F91C4EA84E754F5); // 1.7352737771791877e-2
    const K7: f64 = f64::from_bits(0x3F8C99DE1CA04B32); // 1.396535420284337e-2
    const K8: f64 = f64::from_bits(0x3F87A4D6489AC775); // 1.1544870463101134e-2
    const K9: f64 = f64::from_bits(0x3F842188E828FA97); // 9.829587540591418e-3
    const K10: f64 = f64::from_bits(0x3F803088F95E0E3D); // 7.905073270882851e-3
    const K11: f64 = f64::from_bits(0x3F841E41DCC82CB6); // 9.823336172323099e-3
    const K12: f64 = f64::from_bits(0xBF66F13042EB68EE); // -2.80055452399807e-3
    const K13: f64 = f64::from_bits(0x3F9DBE7F8CBD6964); // 2.9047005620080127e-2
    const K14: f64 = f64::from_bits(0xBFA052C93749A753); // -3.188160705469891e-2
    const K15: f64 = f64::from_bits(0x3FA1EAA418A136A9); // 3.4993293768463225e-2

    // The terms of degree 4 and higher (up to about 2^-8.3 relative to
    // `S(u)`) are evaluated with plain `f64` arithmetic, where the first two
    // steps are in Horner form, so the error of `t` is about 2^-57.5, which
    // is multiplied by u^4 <= 2^-8 in `S(u)`
    let u = u_hi;
    let u2 = u * u;
    let u4 = u2 * u2;
    let t = ((K6 + u * K7) + u2 * (K8 + u * K9))
        + u4 * (((K10 + u * K11) + u2 * (K12 + u * K13)) + u4 * (K14 + u * K15));
    let t = K4 + u * (K5 + u * t);

    // K3 + u * t, where the rounding error of `u * t` (less than 2^-59.8) is
    // multiplied by u^3 <= 2^-6 in `S(u)`, and the next steps use exact
    // products
    let p = u * t;
    let c = F64x2::fast_add11(K3.hi(), p);
    let us = Split::new(u);
    let (h, l) = us.mul_add(c.hi(), c.lo() + K3.lo(), K2);
    let (h, l) = us.mul_add(h, l, K1);
    let (s_hi, s_lo) = us.mul(h, l);

    // S(u_hi + u_lo) ~= S(u_hi) + u_lo * S'(u_hi), where
    // S'(u) = B(u) + u * B'(u), with B(u) = S(u) / u ~= h and
    // B'(u) ~= K2 + 2 * K3 * u + 3 * K4 * u^2 (the next terms are less than
    // 2^-6 relative to S'(u), and u_lo * S'(u) < 2^-58)
    let ds = h + u * (K2.hi() + u * (2.0 * K3.hi() + u * (3.0 * K4)));
    (s_hi, s_lo + u_lo * ds)
}

/// Returns `(hi, lo)` such that `hi + lo ~= atan(x)`, for finite `x`.
///
/// `|lo| <= 2^-14 * |hi|`.
#[inline]
fn atan_core(x: f64) -> (f64, f64) {
    let ax = x.abs();
    let (hi, lo) = if ax <= 1.0 {
        if ax < f64::exp2i_fast(-32) {
            // atan(x) = x - x^3 / 3 + ..., where the second term is less
            // than 2^-65.6 relative to the first one
            (ax, 0.0)
        } else {
            atan_parts(ax, 1.0, ax).finish()
        }
    } else if ax < f64::exp2i_fast(60) {
        // atan(|x|) = π/2 - atan(1 / |x|)
        atan_parts(1.0, ax, 1.0 / ax).finish_with_base(FRAC_PI_2, -1.0)
    } else {
        // atan(|x|) = π/2 - 1 / |x| + ..., where `1 / |x| <= 2^-60` is too
        // small to change the rounded result (π/2, 90 or 0.5 depending on
        // the unit)
        (FRAC_PI_2.hi(), FRAC_PI_2.lo())
    };
    if x.is_sign_negative() {
        (-hi, -lo)
    } else {
        (hi, lo)
    }
}

/// Returns `(hi, lo, edelta)` such that `(hi + lo) * 2^edelta ~= atan2(y, x)`,
/// for finite and non-zero `x` and `y`.
///
/// When `edelta != 0`, `(hi, lo)` is normalized. Otherwise,
/// `|lo| <= 2^-14 * |hi|`.
#[inline]
fn atan2_core(y: f64, x: f64) -> (f64, f64, i32) {
    // atan2(|y|, x) = base + sign * atan(n / d), with n = min(|x|, |y|) and
    // d = max(|x|, |y|):
    // * x > 0, |y| <= |x|: atan2(|y|, x) = atan(n / d)
    // * x > 0, |y| > |x|: atan2(|y|, x) = π/2 - atan(n / d)
    // * x < 0, |y| <= |x|: atan2(|y|, x) = π - atan(n / d)
    // * x < 0, |y| > |x|: atan2(|y|, x) = π/2 + atan(n / d)
    // `(base, sign)` is taken from `BASE_SIGN`, without branches.
    const BASE_SIGN: [(F64x2, f64); 4] = [
        (F64x2::new1(0.0), 1.0),
        (FRAC_PI_2, -1.0),
        (PI, -1.0),
        (FRAC_PI_2, 1.0),
    ];

    let ay = y.abs();
    let ax = x.abs();
    let swap = ay > ax;
    let n = ay.min(ax);
    let d = ay.max(ax);
    let (base, sign) = BASE_SIGN[(usize::from(x.is_sign_negative()) << 1) | usize::from(swap)];

    let q = n / d;
    let (hi, lo, edelta) = if q >= f64::exp2i_fast(-32) {
        // Avoid overflow in `d + c * n` and a subnormal `1 / (d + c * n)`
        // (`n * 0.125` cannot underflow, as `n / d >= 2^-32`)
        let (n, d) = if d.exponent() >= 1021 {
            (n * 0.125, d * 0.125)
        } else {
            (n, d)
        };
        let (hi, lo) = atan_parts(n, d, q).finish_with_base(base, sign);
        (hi, lo, 0)
    } else if base.hi() == 0.0 {
        // atan(n / d) = n / d - (n / d)^3 / 3 + ..., where the second term
        // is less than 2^-65.6 relative to the first one, and the result
        // can be subnormal, so `n / d` is calculated as a double-word scaled
        // by 2^edelta (avoiding subnormals in `F64x2`, which break its
        // accuracy)
        let (n, edelta) = if n.exponent() - d.exponent() < -960 {
            (n * f64::exp2i_fast(106), -106)
        } else {
            (n, 0)
        };
        let z = F64x2::div11(n, d);
        (z.hi(), z.lo(), edelta)
    } else {
        // `n / d < 2^-32` is small compared to the base (π/2 or π), and its
        // rounding error is negligible
        (base.hi(), base.lo() + sign * q, 0)
    };

    if y.is_sign_negative() {
        (-hi, -lo, edelta)
    } else {
        (hi, lo, edelta)
    }
}

/// `atan(n / d)` split as `a + tq + lo` (see `atan_parts`), where `a` is the
/// high part of `atan(c)` (the low part is added to `lo`), and `a >= |tq|` or
/// `a = 0`.
struct AtanParts {
    a: f64,
    tq: f64,
    lo: f64,
}

impl AtanParts {
    /// Returns `(hi, lo)` such that `hi + lo ~= a + tq + lo`.
    #[inline]
    fn finish(self) -> (f64, f64) {
        let s = F64x2::fast_add11(self.a, self.tq);
        (s.hi(), s.lo() + self.lo)
    }

    /// Returns `(hi, lo)` such that `hi + lo ~= base + sign * (a + tq + lo)`,
    /// where `base` is 0, π/2 or π and `sign` is ±1.
    #[inline]
    fn finish_with_base(self, base: F64x2, sign: f64) -> (f64, f64) {
        // |base| >= π/2 > a or base = 0, and |base + sign * a| >= π/4 > |tq|
        // or base + sign * a = ±a, with a >= |tq| or a = 0
        let s1 = F64x2::fast_add11(base.hi(), sign * self.a);
        let s2 = F64x2::fast_add11(s1.hi(), sign * self.tq);
        let lo = (s1.lo() + s2.lo()) + (base.lo() + sign * self.lo);
        (s2.hi(), lo)
    }
}

/// Returns `atan(n / d)` split as `AtanParts`, for `0 <= n <= d < 2^1021`
/// (so `d + c * n < 2^1022` and its reciprocal is normal), given `q ~= n / d`
/// (with a relative error of at most 2^-50).
///
/// The relative error is about 2^-65 (mostly from the polynomial).
#[inline(always)]
fn atan_parts(n: f64, d: f64, q: f64) -> AtanParts {
    // atan(t) ~= t + t^3 * (K2 + K4 * t^2 + K6 * t^4)
    // GENERATE: atan_poly f64 3 1e-30 0.0079
    const K2: f64 = f64::from_bits(0xBFD55555555554FC); // -3.333333333333284e-1
    const K4: f64 = f64::from_bits(0x3FC999999890B46A); // 1.9999999951815833e-1
    const K6: f64 = f64::from_bits(0xBFC248B4D786A7DF); // -1.428438236592333e-1

    // c = i / 64, the closest to n / d
    let (i, c) = atan_index(q);

    // t = (n - c * d) / (d + c * n)
    //
    // `c` has at most 7 significant bits, so its products with the halves of
    // `n` and `d` (with up to 27 bits) are exact.
    // `n - c * d1` is exact: it is `n` when c = 0, and otherwise
    // c * d1 ~= c * d is within a factor of 2 of `n` (Sterbenz lemma).
    let ns = Split::new(n);
    let ds = Split::new(d);
    let num_hi = n - c * ds.hi();
    let num_lo = -(c * ds.lo());
    // `d >= c * n1`
    let den = F64x2::fast_add11(d, c * ns.hi());
    let den_lo = den.lo() + c * ns.lo();

    // t = tq + tl, where `tq` is `t` truncated to 26 bits, so `tq * dh1` and
    // `tq * dh2` are exact (`dh1` and `dh2` are the halves of `den.hi()`).
    // When `num_hi` and `num_lo` nearly cancel, `num_hi - tq * dh1` is not
    // exact, but then `|num_hi| <= 2^-24 * c * d`, so its rounding error is
    // small relative to `atan(c)`.
    let inv = (1.0 / (den.hi() + den_lo)).purify();
    let th = ((num_hi + num_lo) * inv).purify();
    let tq = th.split_hi();
    let dhs = Split::new(den.hi());
    let rem = ((num_hi - tq * dhs.hi()) - tq * dhs.lo()) + (num_lo - tq * den_lo);
    let tl = rem * inv;

    // atan(t) - t
    let t2 = th * th;
    let p = (th * t2) * (K2 + t2 * (K4 + t2 * K6));

    let a = ATAN_TBL[i];
    AtanParts {
        a: a.hi(),
        tq,
        lo: (a.lo() + tl) + p,
    }
}

/// Returns `(i, c)`, where `c = i / 64` is the closest to `q`, for
/// `0 <= q <= 1 + 1/128` (so `i <= 64`, as `ATAN_TBL` requires), and
/// `|q - c| <= 1/128`.
///
/// `c` has at most 7 significant bits.
#[inline]
pub(crate) fn atan_index(q: f64) -> (usize, f64) {
    let (cf, i) = round_i32((q * 64.0).min(64.0));
    (i as usize, cf * (1.0 / 64.0))
}

/// Returns `atan(i / 64)` rounded to `f64`, for `i <= 64`.
#[inline]
pub(crate) fn atan_tbl(i: usize) -> f64 {
    ATAN_TBL[i].hi()
}
