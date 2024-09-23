use super::{print_f32_const, split_hi_lo, FPREC};
use crate::sollya;

pub(crate) fn gen_consts() {
    // π/2
    let mut tmp = rug::Float::with_val(FPREC, rug::float::Constant::Pi) / 2u8;
    let (hi, lo) = split_hi_lo(&mut tmp, 0);
    print_f32_const("FRAC_PI_2_HI", hi);
    print_f32_const("FRAC_PI_2_LO", lo);

    // 3π/4
    let tmp = rug::Float::with_val(FPREC, rug::float::Constant::Pi) * 0.75f32;
    let v = tmp.to_f32();
    print_f32_const("FRAC_3PI_4", v);
}

pub(crate) fn gen_atan_poly() {
    let f = "atan(x) - x";
    let poly_i = [3, 5, 7, 9, 11, 13, 15, 17, 19];
    let range = (-0.001, 1.0);

    sollya::run_and_print_remez_f32(f, range, &poly_i, 0, "K");
}
