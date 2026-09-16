use crate::double::SemiDouble;

impl crate::generic::DivPi for f64 {
    // GENERATE: div_pi::consts f64
    const FRAC_1_PI_EX: SemiDouble<f64> = SemiDouble::with_parts(
        f64::from_bits(0x3FD45F3068000000), // 3.1830988079309464e-1
        f64::from_bits(0x3E3727220A94FE14), // 5.390696036528002e-9
    );
}
