use super::super::{arg_utils, julia, render_const, split_hi_lo, FloatKind};

pub(in super::super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let fkind: FloatKind = arg_utils::parse_1_arg(args)?;
    let aux_prec = fkind.rug_aux_prec();

    let mut out = String::new();

    // 0.5*ln(2π)
    let tmp = (rug::Float::with_val(aux_prec, rug::float::Constant::Pi) * 2u8).ln() / 2u8;
    let (hi, lo) = split_hi_lo(tmp, fkind.float_prec());
    render_const(fkind, "HALF_LN_2_PI_HI", hi, &mut out);
    render_const(fkind, "HALF_LN_2_PI_LO", lo, &mut out);

    Ok(out)
}

pub(in super::super) fn gen_lgamma_poly_1(args: &[&str]) -> Result<String, String> {
    let (fkind, poly_deg, range_end): (_, _, f64) = arg_utils::parse_3_args(args)?;

    let mut out = String::new();

    let func = "SpecialFunctions.lgamma(x + 1)";
    let o = 1.0;
    let range = (0.5 - o, range_end + 0.001 - o);

    julia::run_and_render_remez(fkind, func, range, poly_deg, 0, "K", &[0], &mut out);

    Ok(out)
}

pub(in super::super) fn gen_lgamma_poly_2(args: &[&str]) -> Result<String, String> {
    let (fkind, poly_deg, range_start, range_end): (_, _, f64, f64) =
        arg_utils::parse_4_args(args)?;

    let mut out = String::new();

    let func = "SpecialFunctions.lgamma(x + 2)";
    let o = 2.0;
    let range = (range_start - 0.001 - o, range_end + 0.001 - o);

    julia::run_and_render_remez(fkind, func, range, poly_deg, 0, "K", &[0], &mut out);

    Ok(out)
}

pub(in super::super) fn gen_special_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, poly_deg, range_end): (_, _, f64) = arg_utils::parse_3_args(args)?;

    let mut out = String::new();

    let func = "(SpecialFunctions.gamma(1 / x) / ((1 / x)^(1/x - 0.5) * exp(-1/x) * sqrt(2*BigFloat(pi))) - 1) / x";
    let range = (1.0e-10, 1.0 / (range_end - 0.001));

    julia::run_and_render_remez(fkind, func, range, poly_deg, 0, "K", &[], &mut out);

    Ok(out)
}
