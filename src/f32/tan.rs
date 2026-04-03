impl crate::generic::Tan for f32 {
    #[inline]
    fn tan_poly(x2: Self, x3: Self) -> Self {
        // GENERATE: tan::tan_poly f32 4
        const K3: f32 = f32::from_bits(0x3EAAAA9D); // 3.3333293e-1
        const K5: f32 = f32::from_bits(0x3E088F5D); // 1.3335939e-1
        const K7: f32 = f32::from_bits(0x3D5B0202); // 5.346871e-2
        const K9: f32 = f32::from_bits(0x3CD191B1); // 2.5582166e-2

        horner!(x3, x2, [K3, K5, K7, K9])
    }
}
