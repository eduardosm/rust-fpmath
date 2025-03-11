use super::{print_f32_const, split_hi_lo, FPREC};

pub(crate) fn gen_consts() {
    // 180/π
    let mut tmp = 180u8 / rug::Float::with_val(FPREC, rug::float::Constant::Pi);
    let v = tmp.to_f32();
    let (hi, lo) = split_hi_lo(&mut tmp, 12);
    print_f32_const("RAD_TO_DEG", v);
    print_f32_const("RAD_TO_DEG_HI", hi);
    print_f32_const("RAD_TO_DEG_LO", lo);
}
