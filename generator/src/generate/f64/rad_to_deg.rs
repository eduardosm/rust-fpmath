use super::{print_f64_const, split_hi_lo, FPREC};

pub(crate) fn gen_consts() {
    // 180/π
    let mut tmp = 180u8 / rug::Float::with_val(FPREC, rug::float::Constant::Pi);
    let v = tmp.to_f64();
    let (hi, lo) = split_hi_lo(&mut tmp, 27);
    print_f64_const("RAD_TO_DEG", v);
    print_f64_const("RAD_TO_DEG_HI", hi);
    print_f64_const("RAD_TO_DEG_LO", lo);
}
