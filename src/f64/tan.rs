impl crate::generic::Tan for f64 {
    #[inline]
    fn tan_poly(x2: Self, x3: Self) -> Self {
        // GENERATE: tan::tan_poly f64 9
        const K3: u64 = 0x3FD5555555555575; // 3.333333333333351e-1
        const K5: u64 = 0x3FC111111110D100; // 1.333333333328781e-1
        const K7: u64 = 0x3FABA1BA1BFA9D71; // 5.396825400867556e-2
        const K9: u64 = 0x3F9664F469F9BBAF; // 2.1869486778480936e-2
        const K11: u64 = 0x3F8226E9208FC851; // 8.863278680994709e-3
        const K13: u64 = 0x3F6D6BE6B0CD23F2; // 3.5914903332110388e-3
        const K15: u64 = 0x3F57F26B97A2E986; // 1.4616060930579835e-3
        const K17: u64 = 0x3F425093677176BD; // 5.589217897854142e-4
        const K19: u64 = 0x3F35902BA3288424; // 3.290277992694121e-4

        let k3 = f64::from_bits(K3);
        let k5 = f64::from_bits(K5);
        let k7 = f64::from_bits(K7);
        let k9 = f64::from_bits(K9);
        let k11 = f64::from_bits(K11);
        let k13 = f64::from_bits(K13);
        let k15 = f64::from_bits(K15);
        let k17 = f64::from_bits(K17);
        let k19 = f64::from_bits(K19);

        horner!(x3, x2, [k3, k5, k7, k9, k11, k13, k15, k17, k19])
    }
}
