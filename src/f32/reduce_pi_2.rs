use crate::generic::scalbn_medium;
use crate::traits::{Float as _, FloatConsts as _};

// GENERATE: reduce_pi_2::consts f32
const FRAC_PI_2_HI: u32 = 0x3FC90E00; // 1.5707397e0
const FRAC_PI_2_HIEX: u32 = 0x386D5111; // 5.6580702e-5
const FRAC_PI_2_MI: u32 = 0x386D5000; // 5.657971e-5
const FRAC_PI_2_MIEX: u32 = 0x30885A31; // 9.920936e-10
const FRAC_PI_2_LO: u32 = 0x30885A00; // 9.920882e-10
const FRAC_PI_2_LOEX: u32 = 0x27C234C5; // 5.390303e-15

impl crate::generic::ReducePi2 for f32 {
    #[inline]
    fn frac_pi_2_hi() -> Self {
        f32::from_bits(FRAC_PI_2_HI)
    }

    #[inline]
    fn frac_pi_2_hiex() -> Self {
        f32::from_bits(FRAC_PI_2_HIEX)
    }

    #[inline]
    fn frac_pi_2_mi() -> Self {
        f32::from_bits(FRAC_PI_2_MI)
    }

    #[inline]
    fn frac_pi_2_miex() -> Self {
        f32::from_bits(FRAC_PI_2_MIEX)
    }

    #[inline]
    fn frac_pi_2_lo() -> Self {
        f32::from_bits(FRAC_PI_2_LO)
    }

    #[inline]
    fn frac_pi_2_loex() -> Self {
        f32::from_bits(FRAC_PI_2_LOEX)
    }

    #[inline]
    fn max_reduce_pi_2_medium() -> Self {
        ((1u32 << 9) - 1) as f32 * Self::frac_pi_2()
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
        let fw0 = ((iw as u32) & 0xFFFFFF) as f32;
        let fw1 = (((iw >> 24) as u32) & 0xFFFFFF) as f32 * Self::exp2i_fast(24);
        let fw2 = ((iw >> 48) as u32) as f32 * Self::exp2i_fast(48);

        // add 24-bit chunks into y0, y1
        let mut y0 = ((fw0 + fw1) + fw2).purify();
        let mut y1 = ((fw2 - y0) + fw1) + fw0;

        let scale = i32::from(qe) - 6;
        y0 = scalbn_medium(y0, scale);
        y1 = scalbn_medium(y1, scale);

        if ih == 0 { (y0, y1) } else { (-y0, -y1) }
    }
}
