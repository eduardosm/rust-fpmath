use super::{print_f64_const, split_hi_lo, FPREC};

pub(crate) fn gen_consts() {
    // π/180
    let mut tmp = rug::Float::with_val(FPREC, rug::float::Constant::Pi) / 180u8;
    let v = tmp.to_f64();
    let (hi, lo) = split_hi_lo(&mut tmp, 27);
    print_f64_const("DEG_TO_RAD", v);
    print_f64_const("DEG_TO_RAD_HI", hi);
    print_f64_const("DEG_TO_RAD_LO", lo);
}
