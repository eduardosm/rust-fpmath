use super::super::sollya;
use super::{render_f32_const, split_hi_lo, FPREC};

pub(in super::super) fn gen_consts() -> String {
    let mut out = String::new();

    // sqrt(2)
    let tmp = rug::Float::with_val(FPREC, 2).sqrt();
    let v = tmp.to_f32();
    render_f32_const("SQRT_2", v, &mut out);

    // ln(2)
    let mut tmp = rug::Float::with_val(FPREC, rug::float::Constant::Log2);
    let (hi, lo) = split_hi_lo(&mut tmp, 12);
    render_f32_const("LN_2_HI", hi, &mut out);
    render_f32_const("LN_2_LO", lo, &mut out);

    // 2 / 3
    let mut tmp = rug::Float::with_val(FPREC, 2) / 3u8;
    let (hi, lo) = split_hi_lo(&mut tmp, 0);
    render_f32_const("FRAC_2_3_HI", hi, &mut out);
    render_f32_const("FRAC_2_3_LO", lo, &mut out);

    // 0.4
    let mut tmp = rug::Float::with_val(FPREC, 4) / 10u8;
    let (hi, lo) = split_hi_lo(&mut tmp, 0);
    render_f32_const("FRAC_4_10_HI", hi, &mut out);
    render_f32_const("FRAC_4_10_LO", lo, &mut out);

    out
}

pub(in super::super) fn gen_log_special_poly() -> String {
    let mut out = String::new();

    let f = "(log1p(x) - log(1 - x) - 2 * x) / x";
    let poly_i = [2, 4, 6, 8];
    let range0 = 0.1716;
    let range = (-range0, range0);

    sollya::run_and_render_remez_f32(f, range, &poly_i, 0, "K", &mut out);

    out
}

pub(in super::super) fn gen_log_special_poly_ex() -> String {
    let mut out = String::new();

    let f = "(log1p(x) - log(1 - x) - 2 * x - (2/3) * x^3 - 0.4 * x^5) / x";
    let poly_i = [6, 8, 10];
    let range0 = 0.1716;
    let range = (-range0, range0);

    sollya::run_and_render_remez_f32(f, range, &poly_i, 0, "K", &mut out);

    out
}
