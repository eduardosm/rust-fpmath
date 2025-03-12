use super::super::sollya;
use super::{render_f64_const, split_hi_lo, FPREC};

pub(in super::super) fn gen_consts() -> String {
    let mut out = String::new();

    // π/2
    let mut tmp = rug::Float::with_val(FPREC, rug::float::Constant::Pi) / 2u8;
    let (hi, lo) = split_hi_lo(&mut tmp, 0);
    render_f64_const("FRAC_PI_2_HI", hi, &mut out);
    render_f64_const("FRAC_PI_2_LO", lo, &mut out);

    // 3π/4
    let tmp = rug::Float::with_val(FPREC, rug::float::Constant::Pi) * 0.75f32;
    let v = tmp.to_f64();
    render_f64_const("FRAC_3PI_4", v, &mut out);

    out
}

pub(in super::super) fn gen_atan_poly() -> String {
    let mut out = String::new();

    let f = "atan(x) - x";
    let poly_i = [
        3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41,
    ];
    let range = (-0.001, 1.0);

    sollya::run_and_render_remez_f64(f, range, &poly_i, 0, "K", &mut out);

    out
}
