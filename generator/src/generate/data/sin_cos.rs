use super::super::{arg_utils, render_const, sollya, split_hi_lo, FloatKind};

pub(in super::super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let fkind: FloatKind = arg_utils::parse_1_arg(args)?;
    let aux_prec = fkind.rug_aux_prec();

    let mut out = String::new();

    // 1/6
    let tmp = rug::Float::with_val(aux_prec, 6).recip();
    let (hi, lo) = split_hi_lo(tmp, fkind.split_prec());
    render_const(fkind, "FRAC_1_6_HI", hi, &mut out);
    render_const(fkind, "FRAC_1_6_LO", lo, &mut out);

    Ok(out)
}

pub(in super::super) fn gen_sin_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs) = arg_utils::parse_2_args(args)?;

    let mut out = String::new();

    let func = "sin(x) / x - 1";
    let poly_i = (1..=num_coeffs).map(|i| i * 2).collect::<Vec<_>>();
    let range0 = 0.786; // ~= π/4
    let range = (-range0, range0);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 1, "K", &mut out);

    Ok(out)
}

pub(in super::super) fn gen_sin_poly_ex(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs) = arg_utils::parse_2_args(args)?;

    let mut out = String::new();

    let func = "sin(x) / x - 1 + x^2 / 6";
    let poly_i = (1..=num_coeffs).map(|i| i * 2 + 2).collect::<Vec<_>>();
    let range0 = 0.786; // ~= π/4
    let range = (-range0, range0);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 1, "K", &mut out);

    Ok(out)
}

pub(in super::super) fn gen_cos_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs) = arg_utils::parse_2_args(args)?;

    let mut out = String::new();

    let func = "cos(x) - (1 - 0.5 * x^2)";
    let poly_i = (1..=num_coeffs).map(|i| i * 2 + 2).collect::<Vec<_>>();
    let range0 = 0.786; // ~= π/4
    let range = (-range0, range0);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 0, "K", &mut out);

    Ok(out)
}
