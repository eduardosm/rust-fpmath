use crate::double::SemiDouble;

// GENERATE: sin_cos::consts f32
const FRAC_1_6_HI: u32 = 0x3E2AA000; // 1.6662598e-1
const FRAC_1_6_LO: u32 = 0x382AAAAB; // 4.0690105e-5

impl crate::generic::SinCos for f32 {
    #[inline]
    fn frac_1_6_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(f32::from_bits(FRAC_1_6_HI), f32::from_bits(FRAC_1_6_LO))
    }

    #[inline]
    fn sin_poly(x2: Self, x5: Self) -> (Self, Self) {
        // GENERATE: sin_cos::sin_poly f32 3
        const K3: u32 = 0xBE2AAAA3; // -1.6666655e-1
        const K5: u32 = 0x3C0883AC; // 8.332174e-3
        const K7: u32 = 0xB94CA607; // -1.9516806e-4

        let k3 = f32::from_bits(K3);
        let k5 = f32::from_bits(K5);
        let k7 = f32::from_bits(K7);

        let r = horner!(x5, x2, [k5, k7]);
        (r, k3)
    }

    #[inline]
    fn sin_poly_ex(x2: Self, x5: Self) -> Self {
        // GENERATE: sin_cos::sin_poly_ex f32 2
        const K5: u32 = 0x3C088602; // 8.332731e-3
        const K7: u32 = 0xB94D49A3; // -1.9577755e-4

        let k5 = f32::from_bits(K5);
        let k7 = f32::from_bits(K7);

        horner!(x5, x2, [k5, k7])
    }

    #[inline]
    fn cos_poly(x2: Self, x4: Self) -> Self {
        // GENERATE: sin_cos::cos_poly f32 3
        const K4: u32 = 0x3D2AAAA5; // 4.1666646e-2
        const K6: u32 = 0xBAB60642; // -1.3887363e-3
        const K8: u32 = 0x37CCFFFD; // 2.4437899e-5

        let k4 = f32::from_bits(K4);
        let k6 = f32::from_bits(K6);
        let k8 = f32::from_bits(K8);

        horner!(x4, x2, [k4, k6, k8])
    }
}
