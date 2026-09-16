use super::super::{FloatKind, arg_utils, render_const, sollya};

pub(in super::super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let fkind: FloatKind = arg_utils::parse_1_arg(args)?;
    let sd_fkind = fkind.to_semi_double();

    let mut out = String::new();

    // 1/6
    let tmp = rug::Float::with_val(sd_fkind.rug_aux_prec(), 6).recip();
    render_const(sd_fkind, "FRAC_1_6_EX", tmp, &mut out);

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
