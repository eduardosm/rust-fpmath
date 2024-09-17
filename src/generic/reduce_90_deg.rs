use super::reduce_pi_2::round_fi;
use crate::traits::{CastFrom, CastInto, Float, Int as _, Like, SInt};

pub(crate) trait Reduce90Deg<L = Like<Self>>: Float {
    fn deg_to_rad() -> Self;
    fn deg_to_rad_hi() -> Self;
    fn deg_to_rad_lo() -> Self;

    type SRaw: SInt + CastInto<Self> + CastFrom<Self::Raw>;
}

/// Reduces the angle argument `x` (in degrees) and converts it to
/// radians, returning `(n, y_hi, y_lo)` such as:
/// * `|y_hi| <= π/4`
/// * `|y_lo|` is much smaller than `|y_hi|`
/// * `0 <= n <= 3`
/// * `x = 360*M + 90*n + (y_hi + y_lo)*(180 / π)`
/// * `M` is an integer
pub(crate) fn reduce_90_deg<F: Reduce90Deg>(x: F) -> (u8, F, F) {
    let xabs = x.abs();
    let xexp = x.exponent();
    if xabs <= F::cast_from(45u32) {
        // scale temporarily to avoid subnormal numbers
        let scale = F::exp2i_fast(F::Exp::cast_from(F::MANT_BITS));
        let descale = F::exp2i_fast(-F::Exp::cast_from(F::MANT_BITS));

        let (x_hi, x_lo) = (x * scale).split_hi_lo();
        let y_hi = x_hi * F::deg_to_rad_hi();
        let y_lo = x_hi * F::deg_to_rad_lo() + x_lo * F::deg_to_rad();
        let (y_hi, y_lo) = F::norm_hi_lo_full(y_hi, y_lo);

        (0, y_hi * descale, y_lo * descale)
    } else if xexp <= F::Exp::cast_from(F::MANT_BITS - 4).min(F::Exp::from(31i8)) {
        let (f_n, n) = round_fi(x * (F::one() / F::cast_from(90u32)));
        let ydeg = x - f_n * F::cast_from(90u32);

        let (ydeg_hi, ydeg_lo) = ydeg.split_hi_lo();
        let y_hi = ydeg_hi * F::deg_to_rad_hi();
        let y_lo = ydeg_hi * F::deg_to_rad_lo() + ydeg_lo * F::deg_to_rad();
        let (y_hi, y_lo) = F::norm_hi_lo_full(y_hi, y_lo);

        (n as u8 & 3, y_hi, y_lo)
    } else if xexp < F::Exp::cast_from(F::BITS - 1) {
        let xraw = x.to_raw();

        // ixint = int part of |x|
        // xfrac = frac part of |x|
        let (ixint, xfrac) = if xexp <= F::Exp::cast_from(F::MANT_BITS) {
            let shift: u8 = (F::Exp::cast_from(F::MANT_BITS) - xexp).cast_into();

            let fxint = F::from_raw(xraw & !F::SIGN_MASK & (F::Raw::MAX << shift));
            let ixint = x.mant() >> shift;
            let xfrac = xabs - fxint;

            (ixint, xfrac)
        } else {
            let shift: u8 = (xexp - F::Exp::cast_from(F::MANT_BITS)).cast_into();
            let ixint = x.mant() << shift;

            (ixint, F::ZERO)
        };

        // n = round(trunc(|x|) / 90)
        let n = (ixint + F::Raw::from(45u8)) / F::Raw::from(90u8);
        // rem = trunc(|x|) - n * 90
        let irem = F::SRaw::cast_from(ixint) - F::SRaw::cast_from(F::Raw::from(90u8) * n);
        let frem: F = irem.cast_into();

        let ydeg = xfrac + frem;
        let (ydeg_hi, ydeg_lo) = ydeg.split_hi_lo();
        let y_hi = ydeg_hi * F::deg_to_rad_hi();
        let y_lo = ydeg_hi * F::deg_to_rad_lo() + ydeg_lo * F::deg_to_rad();
        let (y_hi, y_lo) = F::norm_hi_lo_full(y_hi, y_lo);

        let n: u8 = n.cast_into();
        if x.sign() {
            (n.wrapping_neg() & 3, -y_hi, -y_lo)
        } else {
            (n & 3, y_hi, y_lo)
        }
    } else {
        // |x| = xm * 2^xe
        let xm = x.mant();
        let xe = xexp - F::Exp::cast_from(F::MANT_BITS);
        let xe: u16 = xe.cast_into();

        // EXP2_MOD45[i] = mod(2^(i + 3), 45)
        const EXP2_MOD45: [u8; 12] = [1, 2, 4, 8, 16, 32, 19, 38, 31, 17, 34, 23];

        // t = xm * mod(2^xe, 45)
        debug_assert!(xe > 3);
        let t = xm * F::Raw::from(EXP2_MOD45[usize::from((xe - 3) % 12)]);
        // rem45 = mod(|x|, 45)
        let rem45: u16 = (t % F::Raw::from(45u8)).cast_into();
        // rem360 = mod(|x|, 360) = mod(|x|, 45) * 8
        let rem360 = rem45 * 8;

        let mut n = (rem360 / 90) as u8;
        let mut rem90 = (rem360 % 90) as i16;
        if rem90 >= 45 {
            n += 1;
            rem90 -= 90;
        }

        let ydeg = F::cast_from(rem90);

        let y_hi = ydeg * F::deg_to_rad_hi();
        let y_lo = ydeg * F::deg_to_rad_lo();
        let (y_hi, y_lo) = F::norm_hi_lo_full(y_hi, y_lo);

        if x.sign() {
            (n.wrapping_neg() & 3, -y_hi, -y_lo)
        } else {
            (n & 3, y_hi, y_lo)
        }
    }
}
