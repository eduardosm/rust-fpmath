use super::super::sollya;
use super::{render_f64_const, split_hi_lo, FPREC};

pub(in super::super) fn gen_consts() -> String {
    let mut out = String::new();

    // cbrt(2)
    let mut tmp = rug::Float::with_val(FPREC, 2).cbrt();
    let (hi, lo) = split_hi_lo(&mut tmp, 27);
    render_f64_const("CBRT_2_HI", hi, &mut out);
    render_f64_const("CBRT_2_LO", lo, &mut out);

    // cbrt(4)
    let mut tmp = rug::Float::with_val(FPREC, 4).cbrt();
    let (hi, lo) = split_hi_lo(&mut tmp, 27);
    render_f64_const("CBRT_4_HI", hi, &mut out);
    render_f64_const("CBRT_4_LO", lo, &mut out);

    out
}

pub(in super::super) fn gen_inv_cbrt_poly() -> String {
    let mut out = String::new();

    let f = "x^(-1/3)";
    let poly_i = [0, 1, 2, 3, 4, 5];
    let range = (1.0 - 0.001, 2.0 + 0.001);

    sollya::run_and_render_remez_f64(f, range, &poly_i, 0, "K", &mut out);

    out
}
