use super::{F64Like, LikeF64};
use crate::generic::scalbn_medium;

// Generated with `./run-generator.sh f64::reduce_pi_2::consts`
const FRAC_PI_2: u64 = 0x3FF921FB54442D18; // 1.5707963267948966e0
const FRAC_PI_4: u64 = 0x3FE921FB54442D18; // 7.853981633974483e-1
const FRAC_2_PI: u64 = 0x3FE45F306DC9C883; // 6.366197723675814e-1
const FRAC_PI_2_HI: u64 = 0x3FF921FB54400000; // 1.5707963267341256e0
const FRAC_PI_2_HIEX: u64 = 0x3DD0B4611A626331; // 6.077100506506192e-11
const FRAC_PI_2_MI: u64 = 0x3DD0B4611A600000; // 6.077100506303966e-11
const FRAC_PI_2_MIEX: u64 = 0x3BA3198A2E037073; // 2.0222662487959506e-21
const FRAC_PI_2_LO: u64 = 0x3BA3198A2E000000; // 2.0222662487111665e-21
const FRAC_PI_2_LOEX: u64 = 0x397B839A252049C1; // 8.4784276603689e-32

impl<F: F64Like> crate::generic::ReducePi2<LikeF64> for F {
    #[inline]
    fn frac_pi_2() -> Self {
        Self::from_raw(FRAC_PI_2)
    }

    #[inline]
    fn frac_pi_4() -> Self {
        Self::from_raw(FRAC_PI_4)
    }

    #[inline]
    fn frac_2_pi() -> Self {
        Self::from_raw(FRAC_2_PI)
    }

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
        Self::cast_from((1u64 << 20) - 1) * Self::frac_pi_2()
    }

    const REDUCE_PI_2_MEDIUM_TH1: i16 = 16;
    const REDUCE_PI_2_MEDIUM_TH2: i16 = 49;

    type SrcChunks = [u32; 3];

    /// Returns `(x_chunks, e0, jk)`
    fn reduce_pi_2_prepare(x: Self) -> ([u32; 3], i16, usize) {
        let mant = x.mant();
        let x_chunks = [
            ((mant >> 29) as u32),
            (((mant >> 5) & 0x00FF_FFFF) as u32),
            (((mant << 19) & 0x00FF_FFFF) as u32),
        ];
        let e0 = x.exponent() - 23;
        let jk = 4;
        (x_chunks, e0, jk)
    }

    fn reduce_pi_2_compress(qp: &[u64; 20], qe: i16, ih: u32, jz: usize) -> (Self, Self) {
        // iw = sum(qp)
        let mut iw = 0;
        for &qp_i in qp[0..=jz].iter().rev() {
            iw = (iw >> 24) + (u128::from(qp_i) << 48);
        }

        // split iw into 48-bit z
        let fw0 = Self::cast_from((iw as u64) & 0xFFFF_FFFF_FFFF);
        let fw1 = Self::cast_from(((iw >> 48) as u64) & 0xFFFF_FFFF_FFFF) * Self::exp2i_fast(48);
        let fw2 = Self::cast_from((iw >> 96) as u64) * Self::exp2i_fast(96);

        // add 48-bit chunks into y0, y1
        let mut y0 = ((fw0 + fw1) + fw2).purify();
        let mut y1 = ((fw2 - y0) + fw1) + fw0;

        let scale = i32::from(qe) - 48;
        y0 = scalbn_medium(y0, scale);
        y1 = scalbn_medium(y1, scale);

        if ih == 0 {
            (y0, y1)
        } else {
            (-y0, -y1)
        }
    }
}
