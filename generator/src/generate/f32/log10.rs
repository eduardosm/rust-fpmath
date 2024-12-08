use super::{print_f32_const, split_hi_lo, FPREC};

pub(crate) fn gen_consts() {
    // log10(e)
    let mut tmp = rug::Float::with_val(FPREC, 1).exp().log10();
    let (hi, lo) = split_hi_lo(&mut tmp, 12);
    print_f32_const("LOG10_E_HI", hi);
    print_f32_const("LOG10_E_LO", lo);

    // log10(2)
    let mut tmp = rug::Float::with_val(FPREC, 2).log10();
    let (hi, lo) = split_hi_lo(&mut tmp, 12);
    print_f32_const("LOG10_2_HI", hi);
    print_f32_const("LOG10_2_LO", lo);
}
