impl crate::generic::Tan for f64 {
    #[inline]
    fn tan_poly(x2: Self, x3: Self) -> Self {
        // GENERATE: tan::tan_poly f64 9
        const K3: f64 = f64::from_bits(0x3FD5555555555575); // 3.333333333333351e-1
        const K5: f64 = f64::from_bits(0x3FC111111110D100); // 1.333333333328781e-1
        const K7: f64 = f64::from_bits(0x3FABA1BA1BFA9D71); // 5.396825400867556e-2
        const K9: f64 = f64::from_bits(0x3F9664F469F9BBAF); // 2.1869486778480936e-2
        const K11: f64 = f64::from_bits(0x3F8226E9208FC851); // 8.863278680994709e-3
        const K13: f64 = f64::from_bits(0x3F6D6BE6B0CD23F2); // 3.5914903332110388e-3
        const K15: f64 = f64::from_bits(0x3F57F26B97A2E986); // 1.4616060930579835e-3
        const K17: f64 = f64::from_bits(0x3F425093677176BD); // 5.589217897854142e-4
        const K19: f64 = f64::from_bits(0x3F35902BA3288424); // 3.290277992694121e-4

        horner!(x3, x2, [K3, K5, K7, K9, K11, K13, K15, K17, K19])
    }
}
