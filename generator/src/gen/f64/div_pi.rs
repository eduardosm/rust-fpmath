use super::{print_f64_const, split_hi_lo, FPREC};

pub(crate) fn gen_consts() {
    // 1/π
    let mut tmp = 1u8 / rug::Float::with_val(FPREC, rug::float::Constant::Pi);
    let (hi, lo) = split_hi_lo(&mut tmp, 27);
    print_f64_const("FRAC_1_PI_HI", hi);
    print_f64_const("FRAC_1_PI_LO", lo);
}
