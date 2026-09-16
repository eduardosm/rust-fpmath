use crate::double::SemiDouble;

impl crate::generic::Cbrt for f32 {
    // GENERATE: cbrt::consts f32
    const CBRT_2_EX: SemiDouble<f32> = SemiDouble::with_parts(
        f32::from_bits(0x3FA14000), // 1.2597656e0
        f32::from_bits(0x3922F98D), // 1.5542489e-4
    );
    const CBRT_4_EX: SemiDouble<f32> = SemiDouble::with_parts(
        f32::from_bits(0x3FCB2000), // 1.5869141e0
        f32::from_bits(0x39FF529F), // 4.8698948e-4
    );

    #[inline]
    fn exp_mod_3(e: i16) -> i8 {
        (((e + 153) as u16) % 3) as i8
    }

    #[inline]
    fn inv_cbrt_poly(x: Self) -> Self {
        // GENERATE: cbrt::inv_cbrt_poly f32 3
        const K0: f32 = f32::from_bits(0x3FB21939); // 1.3913947e0
        const K1: f32 = f32::from_bits(0xBEF9C752); // -4.8784882e-1
        const K2: f32 = f32::from_bits(0x3DC257A9); // 9.489376e-2

        K0 + horner!(x, x, [K1, K2])
    }
}
