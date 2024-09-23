use super::{print_f32_const, split_hi_lo, FPREC};

pub(crate) fn gen_consts() {
    // log2(10)
    let tmp = rug::Float::with_val(FPREC, 10).log2();
    let v = tmp.to_f32();
    print_f32_const("LOG2_10", v);

    // log10(2)
    let mut tmp = rug::Float::with_val(FPREC, 2).log10();
    let (hi, lo) = split_hi_lo(&mut tmp, 12);
    print_f32_const("LOG10_2_HI", hi);
    print_f32_const("LOG10_2_LO", lo);

    // ln(10)
    let mut tmp = rug::Float::with_val(FPREC, 10).ln();
    let v = tmp.to_f32();
    let (hi, lo) = split_hi_lo(&mut tmp, 12);
    print_f32_const("LN_10", v);
    print_f32_const("LN_10_HI", hi);
    print_f32_const("LN_10_LO", lo);
}
