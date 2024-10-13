use super::{F64Like, LikeF64};
use crate::double::SemiDouble;

// Generated with `./run-generator.sh f64::sin_cos::consts`
const FRAC_1_6_HI: u64 = 0x3FC5555550000000; // 1.666666641831398e-1
const FRAC_1_6_LO: u64 = 0x3E25555555555555; // 2.483526865641276e-9

impl<F: F64Like> crate::generic::SinCos<LikeF64> for F {
    #[inline]
    fn frac_1_6_ex() -> SemiDouble<Self> {
        SemiDouble::with_parts(Self::from_raw(FRAC_1_6_HI), Self::from_raw(FRAC_1_6_LO))
    }

    #[inline]
    fn sin_poly(x2: Self, x5: Self) -> (Self, Self) {
        // Generated with `./run-generator.sh f64::sin_cos::sin_poly`
        const K3: u64 = 0xBFC5555555555549; // -1.6666666666666632e-1
        const K5: u64 = 0x3F8111111110F850; // 8.33333333332234e-3
        const K7: u64 = 0xBF2A01A019C0C181; // -1.98412698297467e-4
        const K9: u64 = 0x3EC71DE3572BCC5E; // 2.7557313669823255e-6
        const K11: u64 = 0xBE5AE5E622FB2468; // -2.505075452524787e-8
        const K13: u64 = 0x3DE5D91CAD15FAA4; // 1.5896580453274288e-10

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
        // Generated with `./run-generator.sh f64::sin_cos::sin_poly_ex`
        const K5: u64 = 0x3F81111111110750; // 8.333333333329002e-3
        const K7: u64 = 0xBF2A01A019D9811A; // -1.9841269834142903e-4
        const K9: u64 = 0x3EC71DE3699EA966; // 2.755731498068118e-6
        const K11: u64 = 0xBE5AE5F2E4324531; // -2.5050935781110632e-8
        const K13: u64 = 0x3DE5DC7074147B84; // 1.590603708120961e-10

        let k5 = Self::from_raw(K5);
        let k7 = Self::from_raw(K7);
        let k9 = Self::from_raw(K9);
        let k11 = Self::from_raw(K11);
        let k13 = Self::from_raw(K13);

        horner!(x5, x2, [k5, k7, k9, k11, k13])
    }

    #[inline]
    fn cos_poly(x2: Self, x4: Self) -> Self {
        // Generated with `./run-generator.sh f64::sin_cos::cos_poly`
        const K4: u64 = 0x3FA555555555554C; // 4.16666666666666e-2
        const K6: u64 = 0xBF56C16C16C15150; // -1.3888888888874025e-3
        const K8: u64 = 0x3EFA01A019CAD16E; // 2.4801587289417634e-5
        const K10: u64 = 0xBE927E4F8066F0E8; // -2.755731433287009e-7
        const K12: u64 = 0x3E21EE9E96F3C559; // 2.0875720523910043e-9
        const K14: u64 = 0xBDA8FAD488D4C37E; // -1.1359500385447443e-11

        let k4 = Self::from_raw(K4);
        let k6 = Self::from_raw(K6);
        let k8 = Self::from_raw(K8);
        let k10 = Self::from_raw(K10);
        let k12 = Self::from_raw(K12);
        let k14 = Self::from_raw(K14);

        horner!(x4, x2, [k4, k6, k8, k10, k12, k14])
    }
}
