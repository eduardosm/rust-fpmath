use super::{F32Like, LikeF32};

impl<F: F32Like> crate::generic::Tan<LikeF32> for F {
    #[inline]
    fn tan_poly(x2: Self, x3: Self) -> Self {
        // GENERATE: tan::tan_poly f32 4
        const K3: u32 = 0x3EAAAA9D; // 3.3333293e-1
        const K5: u32 = 0x3E088F5D; // 1.3335939e-1
        const K7: u32 = 0x3D5B0202; // 5.346871e-2
        const K9: u32 = 0x3CD191B1; // 2.5582166e-2

        let k3 = Self::from_raw(K3);
        let k5 = Self::from_raw(K5);
        let k7 = Self::from_raw(K7);
        let k9 = Self::from_raw(K9);

        horner!(x3, x2, [k3, k5, k7, k9])
    }
}
