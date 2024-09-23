use super::{print_f32_const, split_hi_lo, FPREC};

pub(crate) fn gen_consts() {
    // π
    let mut tmp = rug::Float::with_val(FPREC, rug::float::Constant::Pi);
    let (hi, lo) = split_hi_lo(&mut tmp, 12);
    print_f32_const("PI_HI", hi);
    print_f32_const("PI_LO", lo);
}
