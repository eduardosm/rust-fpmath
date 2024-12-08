use super::{print_f32_const, split_hi_lo, FPREC};

pub(crate) fn gen_consts() {
    // π/2
    let mut tmp = rug::Float::with_val(FPREC, rug::float::Constant::Pi) / 2u8;
    let (hi, hiex) = split_hi_lo(&mut tmp, 9);
    let (mi, miex) = split_hi_lo(&mut tmp, 9);
    let (lo, loex) = split_hi_lo(&mut tmp, 9);
    print_f32_const("FRAC_PI_2_HI", hi);
    print_f32_const("FRAC_PI_2_HIEX", hiex);
    print_f32_const("FRAC_PI_2_MI", mi);
    print_f32_const("FRAC_PI_2_MIEX", miex);
    print_f32_const("FRAC_PI_2_LO", lo);
    print_f32_const("FRAC_PI_2_LOEX", loex);
}