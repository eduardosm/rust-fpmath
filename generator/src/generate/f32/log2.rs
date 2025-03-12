use super::{render_f32_const, split_hi_lo, FPREC};

pub(in super::super) fn gen_consts() -> String {
    let mut out = String::new();

    // log2(e)
    let mut tmp = rug::Float::with_val(FPREC, 1).exp().log2();
    let (hi, lo) = split_hi_lo(&mut tmp, 12);
    render_f32_const("LOG2_E_HI", hi, &mut out);
    render_f32_const("LOG2_E_LO", lo, &mut out);

    out
}
