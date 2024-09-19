use super::reduce_pi_2_large::reduce_pi_2_large;
use crate::traits::{CastInto as _, Float, FloatConsts, Int as _, Like};

pub(crate) trait ReducePi2<L = Like<Self>>: FloatConsts {
    // trunc(π/2)
    fn frac_pi_2_hi() -> Self;
    // π/2 - trunc(π/2)
    fn frac_pi_2_hiex() -> Self;
    // trunc(π/2 - trunc(π/2))
    fn frac_pi_2_mi() -> Self;
    // πpi/2 - trunc(π/2) - trunc(π/2 - trunc(π/2))
    fn frac_pi_2_miex() -> Self;
    // trunc(π/2 - trunc(π/2) - trunc(π/2 - trunc(π/2)))
    fn frac_pi_2_lo() -> Self;
    // π/2 - trunc(π/2 - trunc(π/2) - trunc(π/2 - trunc(π/2)))
    fn frac_pi_2_loex() -> Self;

    fn max_reduce_pi_2_medium() -> Self;

    const REDUCE_PI_2_MEDIUM_TH1: Self::Exp;
    const REDUCE_PI_2_MEDIUM_TH2: Self::Exp;

    type SrcChunks: Sized + AsRef<[u32]>;

    /// Returns `(x_chunks, e0, jk)`
    fn reduce_pi_2_prepare(x: Self) -> (Self::SrcChunks, i16, usize);

    /// Returns `(y0, y1)`
    fn reduce_pi_2_compress(qp: &[u64; 20], qe: i16, ih: u32, jz: usize) -> (Self, Self);
}

/// Rounds `x` and returns it as both float and integer
///
/// `x` must be finite and `0.5 <= abs(x) < 2^min(31, MANT_BITS)`
pub(super) fn round_fi<F: Float>(x: F) -> (F, i32) {
    let e = x.raw_exp();
    if e < F::EXP_OFFSET {
        // 0.5 <= abs(x) < 1
        (F::one().copysign(x), 1 - (i32::from(x.sign()) << 1))
    } else {
        // 1 <= abs(x) < 2^min(31, MANT_BITS)
        let shift = F::RawExp::from(F::MANT_BITS) - (e - F::EXP_OFFSET);
        let imask = F::Raw::MAX << shift;
        let fmask = !imask;
        let xraw = x.to_raw();
        let fpart = xraw & fmask;
        let mut ipart_raw = xraw & !fmask;
        let mut ipart_i: i32 = (x.mant() >> shift).cast_into();
        if fpart > (fmask / F::Raw::TWO) {
            // frac >= 0.5
            ipart_raw += fmask + F::Raw::ONE;
            ipart_i += 1;
        }
        let ipart_f = F::from_raw(ipart_raw);
        if x.sign() {
            ipart_i = -ipart_i;
        }
        (ipart_f, ipart_i)
    }
}

// π/4 < x < MAX_REDUCE_PI_2_MEDIUM
fn reduce_pi_2_medium<F: ReducePi2>(x: F) -> (u8, F, F) {
    // Based on __rem_pio2 (the part after 'medium:') from musl libc
    let (f_n, n) = round_fi(x * F::frac_2_pi());
    // The directed rounding thing from musl has been removed, assume
    // round-to-nearest
    let xexp = x.exponent();
    let mut r = x - f_n * F::frac_pi_2_hi();
    let mut w = f_n * F::frac_pi_2_hiex();
    let mut y0 = (r - w).purify();
    if (xexp - y0.exponent()) > F::REDUCE_PI_2_MEDIUM_TH1 {
        let t = r;
        w = f_n * F::frac_pi_2_mi();
        r = (t - w).purify();
        w = f_n * F::frac_pi_2_miex() - ((t - r) - w);
        y0 = (r - w).purify();
        if (xexp - y0.exponent()) > F::REDUCE_PI_2_MEDIUM_TH2 {
            let t = r;
            w = f_n * F::frac_pi_2_lo();
            r = (t - w).purify();
            w = f_n * F::frac_pi_2_loex() - ((t - r) - w);
            y0 = (r - w).purify();
        }
    }
    let y1 = (r - y0).purify() - w;
    (n as u8 & 3, y0, y1)
}

/// Reduces the angle argument `x`, returning `(n, y_hi, y_lo)`
/// such as:
/// * `|y_hi| <= π/4`
/// * `|y_lo|` is much smaller than `|y_hi|`
/// * `0 <= n <= 3`
/// * `x = 2*π*M + π/2*n + y_hi + y_lo`
/// * `M` is an integer
pub(crate) fn reduce_pi_2<F: ReducePi2>(x: F) -> (u8, F, F) {
    let xabs = x.abs();
    if xabs <= F::frac_pi_4() {
        // reduction not needed
        (0, x, F::ZERO)
    } else if xabs < F::max_reduce_pi_2_medium() {
        reduce_pi_2_medium(x)
    } else {
        let (x_chunks, e0, jk) = F::reduce_pi_2_prepare(x);
        let mut qp: [u64; 20] = [0; 20];
        let (ih, jz, n, qe) = reduce_pi_2_large(x_chunks.as_ref(), e0, jk, &mut qp);
        let (y_hi, y_lo) = F::reduce_pi_2_compress(&qp, qe, ih, jz);

        if x.sign() {
            (n.wrapping_neg() & 3, -y_hi, -y_lo)
        } else {
            (n & 3, y_hi, y_lo)
        }
    }
}
