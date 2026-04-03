use crate::double::NormDouble;

// GENERATE: asin_acos::consts f32
const FRAC_PI_2_HI: f32 = f32::from_bits(0x3FC90FDA); // 1.5707963e0
const FRAC_PI_2_LO: f32 = f32::from_bits(0x33A22169); // 7.54979e-8

impl crate::generic::AsinAcos for f32 {
    #[inline]
    fn frac_pi_2_ex() -> NormDouble<Self> {
        NormDouble::with_parts(FRAC_PI_2_HI, FRAC_PI_2_LO)
    }

    #[inline]
    fn asin_poly(x2: Self) -> Self {
        // GENERATE: asin_acos::asin_poly f32 5
        const K0: f32 = f32::from_bits(0x3E2AAB15); // 1.6666825e-1
        const K2: f32 = f32::from_bits(0x3D99749A); // 7.492943e-2
        const K4: f32 = f32::from_bits(0x3D3B48FF); // 4.572391e-2
        const K6: f32 = f32::from_bits(0x3CBCF147); // 2.3064269e-2
        const K8: f32 = f32::from_bits(0x3D33CAF4); // 4.3894723e-2

        K0 + horner!(x2, x2, [K2, K4, K6, K8])
    }
}
