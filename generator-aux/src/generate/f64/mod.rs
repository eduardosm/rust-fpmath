pub(crate) mod asin_acos;
pub(crate) mod atan;
pub(crate) mod cbrt;
pub(crate) mod div_pi;
pub(crate) mod exp;
pub(crate) mod exp10;
pub(crate) mod exp2;
pub(crate) mod gamma;
pub(crate) mod log;
pub(crate) mod log10;
pub(crate) mod log2;
pub(crate) mod rad_to_deg;
pub(crate) mod reduce_90_deg;
pub(crate) mod reduce_half_mul_pi;
pub(crate) mod reduce_pi_2;
pub(crate) mod sin_cos;
pub(crate) mod tan;

const FPREC: u32 = 1024;

fn split_hi_lo(tmp: &mut rug::Float, n_zeros_in_hi: u8) -> (f64, f64) {
    let hi = tmp.to_f64_round(rug::float::Round::Zero);
    let hi = f64::from_bits(hi.to_bits() & (u64::MAX << n_zeros_in_hi));
    *tmp -= hi;
    let lo = tmp.to_f64();
    (hi, lo)
}

fn print_f64_const<N: std::fmt::Display>(name: N, value: f64) {
    println!(
        "const {name}: u64 = 0x{:016X}; // {value:e}",
        value.to_bits(),
    );
}

pub(crate) fn gen_consts() {
    // π
    let tmp = rug::Float::with_val(FPREC, rug::float::Constant::Pi);
    let v = tmp.to_f64();
    print_f64_const("PI", v);

    // π/2
    let tmp = rug::Float::with_val(FPREC, rug::float::Constant::Pi) / 2u8;
    let v = tmp.to_f64();
    print_f64_const("FRAC_PI_2", v);

    // π/4
    let tmp = rug::Float::with_val(FPREC, rug::float::Constant::Pi) / 4u8;
    let v = tmp.to_f64();
    print_f64_const("FRAC_PI_4", v);

    // 2/π
    let tmp = 2u8 / rug::Float::with_val(FPREC, rug::float::Constant::Pi);
    let v = tmp.to_f64();
    print_f64_const("FRAC_2_PI", v);
}
