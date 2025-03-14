use super::{F64Like, LikeF64};
use crate::double::SemiDouble;

// GENERATE: sin_cos::consts f64
const FRAC_1_6_HI: u64 = 0x3FC5555550000000; // 1.666666641831398e-1
const FRAC_1_6_LO: u64 = 0x3E25555555555555; // 2.483526865641276e-9

impl<F: F64Like> crate::generic::SinCos<LikeF64> for F {
    #[inline]
    fn frac_1_6_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(Self::from_raw(FRAC_1_6_HI), Self::from_raw(FRAC_1_6_LO))
    }

    #[inline]
    fn sin_poly(x2: Self, x5: Self) -> (Self, Self) {
        // GENERATE: sin_cos::sin_poly f64 6
        const K3: u64 = 0xBFC5555555555549; // -1.6666666666666632e-1
        const K5: u64 = 0x3F8111111110F850; // 8.33333333332234e-3
        const K7: u64 = 0xBF2A01A019C0C17F; // -1.9841269829746694e-4
        const K9: u64 = 0x3EC71DE3572BCA6F; // 2.755731366982116e-6
        const K11: u64 = 0xBE5AE5E622F9A50B; // -2.5050754524923147e-8
        const K13: u64 = 0x3DE5D91CACA8DC14; // 1.5896580434788746e-10

        let k3 = Self::from_raw(K3);
        let k5 = Self::from_raw(K5);
        let k7 = Self::from_raw(K7);
        let k9 = Self::from_raw(K9);
        let k11 = Self::from_raw(K11);
        let k13 = Self::from_raw(K13);

        let r = horner!(x5, x2, [k5, k7, k9, k11, k13]);
        (r, k3)
    }

    #[inline]
    fn sin_poly_ex(x2: Self, x5: Self) -> Self {
        // GENERATE: sin_cos::sin_poly_ex f64 5
        const K5: u64 = 0x3F81111111110750; // 8.333333333329002e-3
        const K7: u64 = 0xBF2A01A019D9811B; // -1.9841269834142906e-4
        const K9: u64 = 0x3EC71DE3699EAA4A; // 2.7557314980682144e-6
        const K11: u64 = 0xBE5AE5F2E432F576; // -2.505093578125994e-8
        const K13: u64 = 0x3DE5DC7074471BDE; // 1.5906037089786052e-10

        let k5 = Self::from_raw(K5);
        let k7 = Self::from_raw(K7);
        let k9 = Self::from_raw(K9);
        let k11 = Self::from_raw(K11);
        let k13 = Self::from_raw(K13);

        horner!(x5, x2, [k5, k7, k9, k11, k13])
    }

    #[inline]
    fn cos_poly(x2: Self, x4: Self) -> Self {
        // GENERATE: sin_cos::cos_poly f64 6
        const K4: u64 = 0x3FA555555555554C; // 4.16666666666666e-2
        const K6: u64 = 0xBF56C16C16C15150; // -1.3888888888874025e-3
        const K8: u64 = 0x3EFA01A019CAD16E; // 2.4801587289417634e-5
        const K10: u64 = 0xBE927E4F8066F0EB; // -2.755731433287011e-7
        const K12: u64 = 0x3E21EE9E96F3C746; // 2.0875720523912082e-9
        const K14: u64 = 0xBDA8FAD488D56E3D; // -1.1359500385518062e-11

        let k4 = Self::from_raw(K4);
        let k6 = Self::from_raw(K6);
        let k8 = Self::from_raw(K8);
        let k10 = Self::from_raw(K10);
        let k12 = Self::from_raw(K12);
        let k14 = Self::from_raw(K14);

        horner!(x4, x2, [k4, k6, k8, k10, k12, k14])
    }
}
