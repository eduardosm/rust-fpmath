use super::{F32Like, LikeF32};

impl<F: F32Like> crate::generic::SinhCosh<LikeF32> for F {
    #[inline]
    fn expo2_hi_th() -> Self {
        Self::cast_from(90u32)
    }
}
