use super::{print_f64_const, split_hi_lo, FPREC};
use crate::sollya;

pub(crate) fn gen_consts() {
    // 1/6
    let mut tmp = rug::Float::with_val(FPREC, 6u8).recip();
    let (hi, lo) = split_hi_lo(&mut tmp, 27);
    print_f64_const("FRAC_1_6_HI", hi);
    print_f64_const("FRAC_1_6_LO", lo);
}

pub(crate) fn gen_sin_poly() {
    let f = "sin(x) / x - 1";
    let poly_i = [2, 4, 6, 8, 10, 12];
    let range0 = 0.786; // ~= π/4
    let range = (-range0, range0);

    sollya::run_and_print_remez_f64(f, range, &poly_i, 1, "K");
}

pub(crate) fn gen_sin_poly_ex() {
    let f = "sin(x) / x - 1 + x^2 / 6";
    let poly_i = [4, 6, 8, 10, 12];
    let range0 = 0.786; // ~= π/4
    let range = (-range0, range0);

    sollya::run_and_print_remez_f64(f, range, &poly_i, 1, "K");
}

pub(crate) fn gen_cos_poly() {
    let f = "cos(x) - (1 - 0.5 * x^2)";
    let poly_i = [4, 6, 8, 10, 12, 14];
    let range0 = 0.786; // ~= π/4
    let range = (-range0, range0);

    sollya::run_and_print_remez_f64(f, range, &poly_i, 0, "K");
}
