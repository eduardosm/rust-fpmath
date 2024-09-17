use super::{F64Like, LikeF64};

impl<F: F64Like> crate::generic::SinhCosh<LikeF64> for F {
    fn expo2_hi_th() -> Self {
        Self::cast_from(711u32)
    }
}
