use super::{render_f64_const, split_hi_lo, FPREC};

pub(in super::super) fn gen_consts() -> String {
    let mut out = String::new();

    // π/2
    let mut tmp = rug::Float::with_val(FPREC, rug::float::Constant::Pi) / 2u8;
    let (hi, hiex) = split_hi_lo(&mut tmp, 20);
    let (mi, miex) = split_hi_lo(&mut tmp, 20);
    let (lo, loex) = split_hi_lo(&mut tmp, 20);
    render_f64_const("FRAC_PI_2_HI", hi, &mut out);
    render_f64_const("FRAC_PI_2_HIEX", hiex, &mut out);
    render_f64_const("FRAC_PI_2_MI", mi, &mut out);
    render_f64_const("FRAC_PI_2_MIEX", miex, &mut out);
    render_f64_const("FRAC_PI_2_LO", lo, &mut out);
    render_f64_const("FRAC_PI_2_LOEX", loex, &mut out);

    out
}
