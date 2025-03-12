use super::{render_f32_const, FPREC};

pub(in super::super) fn gen_consts() -> String {
    let mut out = String::new();

    // ln(2)
    let tmp = rug::Float::with_val(FPREC, rug::float::Constant::Log2);
    let v = tmp.to_f32();
    render_f32_const("LN_2", v, &mut out);

    out
}
