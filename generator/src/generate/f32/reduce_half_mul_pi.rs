use super::{render_f32_const, split_hi_lo, FPREC};

pub(in super::super) fn gen_consts() -> String {
    let mut out = String::new();

    // π
    let mut tmp = rug::Float::with_val(FPREC, rug::float::Constant::Pi);
    let (hi, lo) = split_hi_lo(&mut tmp, 12);
    render_f32_const("PI_HI", hi, &mut out);
    render_f32_const("PI_LO", lo, &mut out);

    out
}
