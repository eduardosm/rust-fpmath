use super::super::sollya;
use super::{render_f64_const, split_hi_lo, FPREC};

pub(in super::super) fn gen_consts() -> String {
    let mut out = String::new();

    // log2(e)
    let tmp = rug::Float::with_val(FPREC, 1).exp().log2();
    let v = tmp.to_f64();
    render_f64_const("LOG2_E", v, &mut out);

    // ln(2)
    let mut tmp = rug::Float::with_val(FPREC, 2).ln();
    let (hi, lo) = split_hi_lo(&mut tmp, 27);
    render_f64_const("LN_2_HI", hi, &mut out);
    render_f64_const("LN_2_LO", lo, &mut out);

    out
}

pub(in super::super) fn gen_exp_special_poly() -> String {
    let mut out = String::new();

    let f = "2 - 2 * x / expm1(x) - x";
    let poly_i = [2, 4, 6, 8, 10];
    let range0 = 0.3466; // ~= 0.5*ln(2)
    let range = (-range0, range0);

    sollya::run_and_render_remez_f64(f, range, &poly_i, 0, "K", &mut out);

    out
}

pub(in super::super) fn gen_exp_m1_special_poly() -> String {
    let mut out = String::new();

    let f = "6/x * ((exp(x) + 1)/expm1(x) - 2/x) - 1";
    let poly_i = [2, 4, 6, 8, 10];
    let range0 = 0.3466; // ~= 0.5*ln(2)
    let range = (-range0, range0);

    sollya::run_and_render_remez_f64(f, range, &poly_i, 0, "K", &mut out);

    out
}
