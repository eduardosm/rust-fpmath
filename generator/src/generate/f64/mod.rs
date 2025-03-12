use std::fmt::Write as _;

pub(super) mod asin_acos;
pub(super) mod atan;
pub(super) mod cbrt;
pub(super) mod div_pi;
pub(super) mod exp;
pub(super) mod exp10;
pub(super) mod exp2;
pub(super) mod gamma;
pub(super) mod log;
pub(super) mod log10;
pub(super) mod log2;
pub(super) mod rad_to_deg;
pub(super) mod reduce_90_deg;
pub(super) mod reduce_half_mul_pi;
pub(super) mod reduce_pi_2;
pub(super) mod sin_cos;
pub(super) mod tan;

const FPREC: u32 = 1024;

fn split_hi_lo(tmp: &mut rug::Float, n_zeros_in_hi: u8) -> (f64, f64) {
    let hi = tmp.to_f64_round(rug::float::Round::Zero);
    let hi = f64::from_bits(hi.to_bits() & (u64::MAX << n_zeros_in_hi));
    *tmp -= hi;
    let lo = tmp.to_f64();
    (hi, lo)
}

fn render_f64_const<N: std::fmt::Display>(name: N, value: f64, out: &mut String) {
    writeln!(
        out,
        "const {name}: u64 = 0x{:016X}; // {value:e}",
        value.to_bits(),
    )
    .unwrap();
}

pub(super) fn gen_consts() -> String {
    let mut out = String::new();

    // π
    let tmp = rug::Float::with_val(FPREC, rug::float::Constant::Pi);
    let v = tmp.to_f64();
    render_f64_const("PI", v, &mut out);

    // π/2
    let tmp = rug::Float::with_val(FPREC, rug::float::Constant::Pi) / 2u8;
    let v = tmp.to_f64();
    render_f64_const("FRAC_PI_2", v, &mut out);

    // π/4
    let tmp = rug::Float::with_val(FPREC, rug::float::Constant::Pi) / 4u8;
    let v = tmp.to_f64();
    render_f64_const("FRAC_PI_4", v, &mut out);

    // 2/π
    let tmp = 2u8 / rug::Float::with_val(FPREC, rug::float::Constant::Pi);
    let v = tmp.to_f64();
    render_f64_const("FRAC_2_PI", v, &mut out);

    out
}
