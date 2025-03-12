use super::super::julia;
use super::{render_f32_const, split_hi_lo, FPREC};

pub(in super::super) fn gen_consts() -> String {
    let mut out = String::new();

    // 0.5*ln(2π)
    let mut tmp = (rug::Float::with_val(FPREC, rug::float::Constant::Pi) * 2u8).ln() / 2u8;
    let (hi, lo) = split_hi_lo(&mut tmp, 0);
    render_f32_const("HALF_LN_2_PI_HI", hi, &mut out);
    render_f32_const("HALF_LN_2_PI_LO", lo, &mut out);

    out
}

pub(in super::super) fn gen_lgamma_poly_1() -> String {
    let mut out = String::new();

    let f = "SpecialFunctions.lgamma(x + 1)";
    let poly_deg = 12;
    let o = 1.0;
    let range = (0.5 - o, 1.201 - o);

    julia::run_and_render_remez_f32(f, range, poly_deg, 0, "K", &[0], &mut out);

    out
}

pub(in super::super) fn gen_lgamma_poly_2() -> String {
    let mut out = String::new();

    let f = "SpecialFunctions.lgamma(x + 2)";
    let poly_deg = 12;
    let o = 2.0;
    let range = (1.199 - o, 2.301 - o);

    julia::run_and_render_remez_f32(f, range, poly_deg, 0, "K", &[0], &mut out);

    out
}

pub(in super::super) fn gen_special_poly() -> String {
    let mut out = String::new();

    let f = "(SpecialFunctions.gamma(1 / x) / ((1 / x)^(1/x - 0.5) * exp(-1/x) * sqrt(2*BigFloat(pi))) - 1) / x";
    let poly_deg = 8;
    let range = (1.0e-10, 1.0 / 2.299);

    julia::run_and_render_remez_f32(f, range, poly_deg, 0, "K", &[], &mut out);

    out
}
