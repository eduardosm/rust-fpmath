use crate::traits::{CastFrom as _, CastInto as _, FloatConsts, Int as _, Like};

pub(crate) trait ReduceHalfMulPi<L = Like<Self>>: FloatConsts {
    fn pi_hi() -> Self;
    fn pi_lo() -> Self;
}

/// Reduces the angle argument `x` (in radians/π) and converts it to
///  radians, returning `(n, y_hi, y_lo)` such as:
/// * `|y_hi| <= π/4`
/// * `|y_lo|` is much smaller than `|y_hi|`
/// * `0 <= n <= 3`
/// * `x = 2*M + 0.5*n + (y_hi + y_lo)/π`
/// * `M` is an integer
pub(crate) fn reduce_half_mul_pi<F: ReduceHalfMulPi>(x: F) -> (u8, F, F) {
    let xexp = x.exponent();
    if xexp < -F::Exp::TWO {
        // |x| < 0.25

        // scale temporarily to avoid temporary subnormal numbers
        let scale = F::exp2i_fast(F::Exp::cast_from(F::MANT_BITS));
        let descale = F::exp2i_fast(-F::Exp::cast_from(F::MANT_BITS));

        let (x_hi, x_lo) = (x * scale).split_hi_lo();
        let y_hi = x_hi * F::pi_hi();
        let y_lo = x_hi * F::pi_lo() + x_lo * F::pi();
        let (y_hi, y_lo) = F::norm_hi_lo_full(y_hi, y_lo);

        (0, y_hi * descale, y_lo * descale)
    } else if xexp == -F::Exp::TWO {
        // 0.25 <= abs(x) < 0.5
        let fpart = x - F::half().copysign(x);

        let n = if x.sign() { 3 } else { 1 };

        let (fpart_hi, fpart_lo) = fpart.split_hi_lo();
        let y_hi = fpart_hi * F::pi_hi();
        let y_lo = fpart_hi * F::pi_lo() + fpart_lo * F::pi();
        let (y_hi, y_lo) = F::norm_hi_lo_full(y_hi, y_lo);

        (n, y_hi, y_lo)
    } else if xexp < F::Exp::cast_from(F::MANT_BITS) {
        // Split x = fpart + ipart * 0.5
        let shift: u8 = (F::Exp::cast_from(F::MANT_BITS) - (xexp + F::Exp::ONE)).cast_into();
        let fmask = !(F::Raw::MAX << shift);
        let xraw = x.to_raw();
        let fpart = xraw & fmask;
        let ipart = xraw & !fmask;
        let mut ipart_i = x.mant() >> shift;
        let mut fpart_f = x - F::from_raw(ipart);
        // Round to nearest
        if fpart > (fmask / F::Raw::TWO) {
            fpart_f = fpart_f - F::half().copysign(x);
            ipart_i += F::Raw::ONE;
        }

        let mut n: u8 = ipart_i.cast_into();
        if x.sign() {
            n = n.wrapping_neg();
        }

        let (fpart_hi, fpart_lo) = fpart_f.split_hi_lo();
        let y_hi = fpart_hi * F::pi_hi();
        let y_lo = fpart_hi * F::pi_lo() + fpart_lo * F::pi();
        let (y_hi, y_lo) = F::norm_hi_lo_full(y_hi, y_lo);

        (n & 3, y_hi, y_lo)
    } else if xexp == F::Exp::cast_from(F::MANT_BITS) {
        // The lowest bit of the integer part is zero
        let n: u8 = ((x.to_raw() & F::Raw::ONE) << 1).cast_into();
        if x.sign() {
            (n, F::ZERO, F::ZERO)
        } else {
            (n.wrapping_neg() & 3, F::ZERO, F::ZERO)
        }
    } else {
        // The two lowest bits of the integer part are zero
        (0, F::ZERO, F::ZERO)
    }
}
