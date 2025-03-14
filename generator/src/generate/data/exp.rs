use super::super::{arg_utils, render_const, sollya, split_hi_lo, FloatKind};

pub(in super::super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let fkind: FloatKind = arg_utils::parse_1_arg(args)?;
    let aux_prec = fkind.rug_aux_prec();

    let mut out = String::new();

    // log2(e)
    let tmp = rug::Float::with_val(aux_prec, 1).exp().log2();
    render_const(fkind, "LOG2_E", tmp, &mut out);

    // ln(2)
    let tmp = rug::Float::with_val(aux_prec, 2).ln();
    let (hi, lo) = split_hi_lo(tmp, fkind.split_prec());
    render_const(fkind, "LN_2_HI", hi, &mut out);
    render_const(fkind, "LN_2_LO", lo, &mut out);

    Ok(out)
}

pub(in super::super) fn gen_exp_special_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs) = arg_utils::parse_2_args(args)?;

    let mut out = String::new();

    let func = "2 - 2 * x / expm1(x) - x";
    let poly_i = (1..=num_coeffs).map(|i| i * 2).collect::<Vec<_>>();
    let range0 = 0.3466; // ~= 0.5*ln(2)
    let range = (-range0, range0);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 0, "K", &mut out);

    Ok(out)
}

pub(in super::super) fn gen_exp_m1_special_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs) = arg_utils::parse_2_args(args)?;

    let mut out = String::new();

    let func = "6/x * ((exp(x) + 1)/expm1(x) - 2/x) - 1";
    let poly_i = (1..=num_coeffs).map(|i| i * 2).collect::<Vec<_>>();
    let range0 = 0.3466; // ~= 0.5*ln(2)
    let range = (-range0, range0);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 0, "K", &mut out);

    Ok(out)
}
