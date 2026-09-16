use crate::double::SemiDouble;

impl crate::generic::Log10 for f64 {
    // GENERATE: log10::consts f64
    const LOG10_E_EX: SemiDouble<f64> = SemiDouble::with_parts(
        f64::from_bits(0x3FDBCB7B10000000),
        f64::from_bits(0x3E349B9438CA9AAE),
    ); // 4.3429448190325182765113e-1
    const LOG10_2_HI: f64 = f64::from_bits(0x3FD3441350000000); // 3.010299950838089e-1
    const LOG10_2_LO: f64 = f64::from_bits(0x3E03EF3FDE623E25); // 5.801722962879576e-10
}
