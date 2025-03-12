use super::{render_f64_const, split_hi_lo, FPREC};

pub(in super::super) fn gen_consts() -> String {
    let mut out = String::new();

    // 180/π
    let mut tmp = 180u8 / rug::Float::with_val(FPREC, rug::float::Constant::Pi);
    let v = tmp.to_f64();
    let (hi, lo) = split_hi_lo(&mut tmp, 27);
    render_f64_const("RAD_TO_DEG", v, &mut out);
    render_f64_const("RAD_TO_DEG_HI", hi, &mut out);
    render_f64_const("RAD_TO_DEG_LO", lo, &mut out);

    out
}
