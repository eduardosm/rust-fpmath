use super::{render_f64_const, split_hi_lo, FPREC};

pub(in super::super) fn gen_consts() -> String {
    let mut out = String::new();

    // log2(10)
    let tmp = rug::Float::with_val(FPREC, 10).log2();
    let v = tmp.to_f64();
    render_f64_const("LOG2_10", v, &mut out);

    // log10(2)
    let mut tmp = rug::Float::with_val(FPREC, 2).log10();
    let (hi, lo) = split_hi_lo(&mut tmp, 27);
    render_f64_const("LOG10_2_HI", hi, &mut out);
    render_f64_const("LOG10_2_LO", lo, &mut out);

    // ln(10)
    let mut tmp = rug::Float::with_val(FPREC, 10).ln();
    let v = tmp.to_f64();
    let (hi, lo) = split_hi_lo(&mut tmp, 27);
    render_f64_const("LN_10", v, &mut out);
    render_f64_const("LN_10_HI", hi, &mut out);
    render_f64_const("LN_10_LO", lo, &mut out);

    out
}
