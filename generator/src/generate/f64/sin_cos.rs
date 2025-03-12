use super::super::sollya;
use super::{render_f64_const, split_hi_lo, FPREC};

pub(in super::super) fn gen_consts() -> String {
    let mut out = String::new();

    // 1/6
    let mut tmp = rug::Float::with_val(FPREC, 6u8).recip();
    let (hi, lo) = split_hi_lo(&mut tmp, 27);
    render_f64_const("FRAC_1_6_HI", hi, &mut out);
    render_f64_const("FRAC_1_6_LO", lo, &mut out);

    out
}

pub(in super::super) fn gen_sin_poly() -> String {
    let mut out = String::new();

    let f = "sin(x) / x - 1";
    let poly_i = [2, 4, 6, 8, 10, 12];
    let range0 = 0.786; // ~= π/4
    let range = (-range0, range0);

    sollya::run_and_render_remez_f64(f, range, &poly_i, 1, "K", &mut out);

    out
}

pub(in super::super) fn gen_sin_poly_ex() -> String {
    let mut out = String::new();

    let f = "sin(x) / x - 1 + x^2 / 6";
    let poly_i = [4, 6, 8, 10, 12];
    let range0 = 0.786; // ~= π/4
    let range = (-range0, range0);

    sollya::run_and_render_remez_f64(f, range, &poly_i, 1, "K", &mut out);

    out
}

pub(in super::super) fn gen_cos_poly() -> String {
    let mut out = String::new();

    let f = "cos(x) - (1 - 0.5 * x^2)";
    let poly_i = [4, 6, 8, 10, 12, 14];
    let range0 = 0.786; // ~= π/4
    let range = (-range0, range0);

    sollya::run_and_render_remez_f64(f, range, &poly_i, 0, "K", &mut out);

    out
}
