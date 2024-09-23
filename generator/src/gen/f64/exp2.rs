use super::{print_f64_const, split_hi_lo, FPREC};

pub(crate) fn gen_consts() {
    // ln(2)
    let mut tmp = rug::Float::with_val(FPREC, rug::float::Constant::Log2);
    let v = tmp.to_f64();
    let (hi, lo) = split_hi_lo(&mut tmp, 27);
    print_f64_const("LN_2", v);
    print_f64_const("LN_2_HI", hi);
    print_f64_const("LN_2_LO", lo);
}
