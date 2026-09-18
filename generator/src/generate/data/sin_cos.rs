use super::super::{arg_utils, sollya};

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
