use crate::double::{NormDouble, SemiDouble};
use crate::traits::{CastFrom as _, CastInto as _, FloatConsts, Int as _, Like};

pub(crate) trait ReduceHalfMulPi<L = Like<Self>>: FloatConsts {
    fn pi_ex() -> SemiDouble<Self>;
}

/// Reduces the angle argument `x` (in radians/π) and converts it to
///  radians, returning `(n, y_hi, y_lo)` such as:
/// * `|y_hi| <= π/4`
/// * `|y_lo|` is much smaller than `|y_hi|`
/// * `0 <= n <= 3`
/// * `x = 2*M + 0.5*n + (y_hi + y_lo)/π`
/// * `M` is an integer
pub(crate) fn reduce_half_mul_pi<F: ReduceHalfMulPi>(x: F) -> (u8, NormDouble<F>) {
    let xexp = x.exponent();
    if xexp < -F::Exp::TWO {
        // |x| < 0.25

        // scale temporarily to avoid temporary subnormal numbers
        let scale = F::exp2i_fast(F::Exp::cast_from(F::MANT_BITS));
        let descale = F::exp2i_fast(-F::Exp::cast_from(F::MANT_BITS));

        let sx = SemiDouble::new(x * scale);
        let y = (sx * F::pi_ex()).pmul1(descale).to_norm();

        (0, y)
    } else if xexp == -F::Exp::TWO {
        // 0.25 <= abs(x) < 0.5
        let fpart = x - F::half().copysign(x);

        let n = if x.sign() { 3 } else { 1 };

        let fpart = SemiDouble::new(fpart);
        let y = (fpart * F::pi_ex()).to_norm();

        (n, y)
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

        let fpart = SemiDouble::new(fpart_f);
        let y = (fpart * F::pi_ex()).to_norm();

        (n & 3, y)
    } else if xexp == F::Exp::cast_from(F::MANT_BITS) {
        // The lowest bit of the integer part is zero
        let n: u8 = ((x.to_raw() & F::Raw::ONE) << 1).cast_into();
        if x.sign() {
            (n, NormDouble::ZERO)
        } else {
            (n.wrapping_neg() & 3, NormDouble::ZERO)
        }
    } else {
        // The two lowest bits of the integer part are zero
        (0, NormDouble::ZERO)
    }
}
