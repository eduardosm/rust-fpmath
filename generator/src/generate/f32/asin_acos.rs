use super::super::sollya;
use super::{render_f32_const, split_hi_lo, FPREC};

pub(in super::super) fn gen_consts() -> String {
    let mut out = String::new();

    // π/2
    let mut tmp = rug::Float::with_val(FPREC, rug::float::Constant::Pi) / 2u8;
    let (hi, lo) = split_hi_lo(&mut tmp, 0);
    render_f32_const("FRAC_PI_2_HI", hi, &mut out);
    render_f32_const("FRAC_PI_2_LO", lo, &mut out);

    out
}

pub(in super::super) fn gen_asin_poly() -> String {
    let mut out = String::new();

    let f = "asin(x) - x";
    let poly_i = [3, 5, 7, 9, 11];
    let range = (-0.001, 0.501);

    sollya::run_and_render_remez_f32(f, range, &poly_i, -3, "K", &mut out);

    out
}
