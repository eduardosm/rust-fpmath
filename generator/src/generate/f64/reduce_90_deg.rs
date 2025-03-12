use super::{render_f64_const, split_hi_lo, FPREC};

pub(in super::super) fn gen_consts() -> String {
    let mut out = String::new();

    // π/180
    let mut tmp = rug::Float::with_val(FPREC, rug::float::Constant::Pi) / 180u8;
    let v = tmp.to_f64();
    let (hi, lo) = split_hi_lo(&mut tmp, 27);
    render_f64_const("DEG_TO_RAD", v, &mut out);
    render_f64_const("DEG_TO_RAD_HI", hi, &mut out);
    render_f64_const("DEG_TO_RAD_LO", lo, &mut out);

    out
}
