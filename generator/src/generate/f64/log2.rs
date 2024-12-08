use super::{print_f64_const, split_hi_lo, FPREC};

pub(crate) fn gen_consts() {
    // log2(e)
    let mut tmp = rug::Float::with_val(FPREC, 1).exp().log2();
    let (hi, lo) = split_hi_lo(&mut tmp, 27);
    print_f64_const("LOG2_E_HI", hi);
    print_f64_const("LOG2_E_LO", lo);
}
