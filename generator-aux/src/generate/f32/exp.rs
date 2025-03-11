use super::{print_f32_const, split_hi_lo, FPREC};
use crate::sollya;

pub(crate) fn gen_consts() {
    // log2(e)
    let tmp = rug::Float::with_val(FPREC, 1).exp().log2();
    let v = tmp.to_f32();
    print_f32_const("LOG2_E", v);

    // ln(2)
    let mut tmp = rug::Float::with_val(FPREC, 2).ln();
    let (hi, lo) = split_hi_lo(&mut tmp, 12);
    print_f32_const("LN_2_HI", hi);
    print_f32_const("LN_2_LO", lo);
}

pub(crate) fn gen_exp_special_poly() {
    let f = "2 - 2 * x / expm1(x) - x";
    let poly_i = [2, 4];
    let range0 = 0.3466; // ~= 0.5*ln(2)
    let range = (-range0, range0);

    sollya::run_and_print_remez_f32(f, range, &poly_i, 0, "K");
}

pub(crate) fn gen_exp_m1_special_poly() {
    let f = "6/x * ((exp(x) + 1)/expm1(x) - 2/x) - 1";
    let poly_i = [2, 4];
    let range0 = 0.3466; // ~= 0.5*ln(2)
    let range = (-range0, range0);

    sollya::run_and_print_remez_f32(f, range, &poly_i, 0, "K");
}
