use super::{render_f32_const, split_hi_lo, FPREC};

pub(in super::super) fn gen_consts() -> String {
    let mut out = String::new();

    // 180/π
    let mut tmp = 180u8 / rug::Float::with_val(FPREC, rug::float::Constant::Pi);
    let v = tmp.to_f32();
    let (hi, lo) = split_hi_lo(&mut tmp, 12);
    render_f32_const("RAD_TO_DEG", v, &mut out);
    render_f32_const("RAD_TO_DEG_HI", hi, &mut out);
    render_f32_const("RAD_TO_DEG_LO", lo, &mut out);

    out
}
