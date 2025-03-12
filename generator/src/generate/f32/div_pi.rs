use super::{render_f32_const, split_hi_lo, FPREC};

pub(in super::super) fn gen_consts() -> String {
    let mut out = String::new();

    // 1/π
    let mut tmp = 1u8 / rug::Float::with_val(FPREC, rug::float::Constant::Pi);
    let (hi, lo) = split_hi_lo(&mut tmp, 12);
    render_f32_const("FRAC_1_PI_HI", hi, &mut out);
    render_f32_const("FRAC_1_PI_LO", lo, &mut out);

    out
}
