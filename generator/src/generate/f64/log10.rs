use super::{render_f64_const, split_hi_lo, FPREC};

pub(in super::super) fn gen_consts() -> String {
    let mut out = String::new();

    // log10(e)
    let mut tmp = rug::Float::with_val(FPREC, 1).exp().log10();
    let (hi, lo) = split_hi_lo(&mut tmp, 27);
    render_f64_const("LOG10_E_HI", hi, &mut out);
    render_f64_const("LOG10_E_LO", lo, &mut out);

    // log10(2)
    let mut tmp = rug::Float::with_val(FPREC, 2).log10();
    let (hi, lo) = split_hi_lo(&mut tmp, 27);
    render_f64_const("LOG10_2_HI", hi, &mut out);
    render_f64_const("LOG10_2_LO", lo, &mut out);

    out
}
