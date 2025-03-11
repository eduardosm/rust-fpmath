use super::{F64Like, LikeF64};

impl<F: F64Like> crate::generic::Tan<LikeF64> for F {
    #[inline]
    fn tan_poly(x2: Self, x3: Self) -> Self {
        // GENERATE: other f64::tan::tan_poly
        const K3: u64 = 0x3FD5555555555575; // 3.333333333333351e-1
        const K5: u64 = 0x3FC111111110D100; // 1.333333333328781e-1
        const K7: u64 = 0x3FABA1BA1BFA9D71; // 5.396825400867556e-2
        const K9: u64 = 0x3F9664F469F9BBAF; // 2.1869486778480936e-2
        const K11: u64 = 0x3F8226E9208FC851; // 8.863278680994709e-3
        const K13: u64 = 0x3F6D6BE6B0CD23F2; // 3.5914903332110388e-3
        const K15: u64 = 0x3F57F26B97A2E986; // 1.4616060930579835e-3
        const K17: u64 = 0x3F425093677176BD; // 5.589217897854142e-4
        const K19: u64 = 0x3F35902BA3288424; // 3.290277992694121e-4

        let k3 = Self::from_raw(K3);
        let k5 = Self::from_raw(K5);
        let k7 = Self::from_raw(K7);
        let k9 = Self::from_raw(K9);
        let k11 = Self::from_raw(K11);
        let k13 = Self::from_raw(K13);
        let k15 = Self::from_raw(K15);
        let k17 = Self::from_raw(K17);
        let k19 = Self::from_raw(K19);

        horner!(x3, x2, [k3, k5, k7, k9, k11, k13, k15, k17, k19])
    }
}
