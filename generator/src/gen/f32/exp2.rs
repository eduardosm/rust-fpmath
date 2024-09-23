use super::{print_f32_const, FPREC};

pub(crate) fn gen_consts() {
    // ln(2)
    let tmp = rug::Float::with_val(FPREC, rug::float::Constant::Log2);
    let v = tmp.to_f32();
    print_f32_const("LN_2", v);
}
