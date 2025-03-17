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

pub(in super::super) fn gen_lgamma_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, poly_deg, offset, range_start, range_end): (_, i32, u8, f64, f64) =
        arg_utils::parse_5_args(args)?;

    if !matches!(offset, 1 | 2) {
        return Err("offset must be 1 or 2".into());
    }
    let offset = f64::from(offset);

    let mut out = String::new();

    let func = format!("SpecialFunctions.lgamma(x + {offset}) / x");
    let wfunc = "1 / fx";
    let range = (range_start - offset, range_end - offset);

    julia::run_and_render_remez(fkind, &func, wfunc, range, poly_deg - 1, 1, "K", &mut out);

    Ok(out)
}

pub(in super::super) fn gen_special_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, poly_deg, range_end): (_, _, f64) = arg_utils::parse_3_args(args)?;

    let mut out = String::new();

    // (SpecialFunctions.gamma(1 / x) / ((1 / x)^(1/x - 0.5) * exp(-1/x) * sqrt(2*BigFloat(pi))) - 1) / x
    let func = "(exp(SpecialFunctions.lgamma(1 / x) + (1/x - 0.5) * log(x) + 1 / x - log(sqrt(2*BigFloat(pi)))) - 1) / x";
    let wfunc = "1";
    let range = (1.0e-100, 1.0 / (range_end - 0.0001));

    julia::run_and_render_remez(fkind, func, wfunc, range, poly_deg, 0, "K", &mut out);

    Ok(out)
}
