use super::{F32Like, LikeF32};
use crate::generic::scalbn_medium;
use crate::traits::FloatConsts;

// Generated with `./run-generator.sh f32::reduce_pi_2::consts`
const FRAC_PI_2_HI: u32 = 0x3FC90E00; // 1.5707397e0
const FRAC_PI_2_HIEX: u32 = 0x386D5111; // 5.6580702e-5
const FRAC_PI_2_MI: u32 = 0x386D5000; // 5.657971e-5
const FRAC_PI_2_MIEX: u32 = 0x30885A31; // 9.920936e-10
const FRAC_PI_2_LO: u32 = 0x30885A00; // 9.920882e-10
const FRAC_PI_2_LOEX: u32 = 0x27C234C5; // 5.390303e-15

impl<F: F32Like + FloatConsts> crate::generic::ReducePi2<LikeF32> for F {
    #[inline]
    fn frac_pi_2_hi() -> Self {
        Self::from_raw(FRAC_PI_2_HI)
    }

    #[inline]
    fn frac_pi_2_hiex() -> Self {
        Self::from_raw(FRAC_PI_2_HIEX)
    }

    #[inline]
    fn frac_pi_2_mi() -> Self {
        Self::from_raw(FRAC_PI_2_MI)
    }

    #[inline]
    fn frac_pi_2_miex() -> Self {
        Self::from_raw(FRAC_PI_2_MIEX)
    }

    #[inline]
    fn frac_pi_2_lo() -> Self {
        Self::from_raw(FRAC_PI_2_LO)
    }

    #[inline]
    fn frac_pi_2_loex() -> Self {
        Self::from_raw(FRAC_PI_2_LOEX)
    }

    fn max_reduce_pi_2_medium() -> Self {
        F::cast_from((1u32 << 9) - 1) * Self::frac_pi_2()
    }

    const REDUCE_PI_2_MEDIUM_TH1: i16 = 8;
    const REDUCE_PI_2_MEDIUM_TH2: i16 = 20;

    type SrcChunks = [u32; 1];

    /// Returns `(x_chunks, e0, jk)`
    fn reduce_pi_2_prepare(x: Self) -> ([u32; 1], i16, usize) {
        let mant = x.mant();
        let x_chunks = [mant];
        let e0 = x.exponent() - 23;
        let jk = 3;
        (x_chunks, e0, jk)
    }

    fn reduce_pi_2_compress(qp: &[u64], qe: i16, ih: u32) -> (Self, Self) {
        // iw = sum(qp)
        let mut iw = 0;
        for &qp_i in qp.iter().rev() {
            iw = (iw >> 24) + (qp_i << 6);
        }

        // split iw into 24-bit chunks
        let fw0 = Self::cast_from((iw as u32) & 0xFFFFFF);
        let fw1 = Self::cast_from(((iw >> 24) as u32) & 0xFFFFFF) * Self::exp2i_fast(24);
        let fw2 = Self::cast_from((iw >> 48) as u32) * Self::exp2i_fast(48);

        // add 24-bit chunks into y0, y1
        let mut y0 = ((fw0 + fw1) + fw2).purify();
        let mut y1 = ((fw2 - y0) + fw1) + fw0;

        let scale = i32::from(qe) - 6;
        y0 = scalbn_medium(y0, scale);
        y1 = scalbn_medium(y1, scale);

        if ih == 0 {
            (y0, y1)
        } else {
            (-y0, -y1)
        }
    }
}
