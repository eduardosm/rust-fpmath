use super::super::{FloatKind, arg_utils, render_const, sollya, split_hi_lo};

pub(in super::super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let fkind: FloatKind = arg_utils::parse_1_arg(args)?;
    let nd_fkind = fkind.to_double();

    let mut out = String::new();

    // sqrt(2)
    let tmp = rug::Float::with_val(fkind.rug_aux_prec(), 2).sqrt();
    render_const(fkind, "SQRT_2", tmp, &mut out);

    // ln(2)
    let tmp = rug::Float::with_val(fkind.rug_aux_prec(), rug::float::Constant::Log2);
    let (hi, lo) = split_hi_lo(tmp, fkind.split_prec());
    render_const(fkind, "LN_2_HI", hi, &mut out);
    render_const(fkind, "LN_2_LO", lo, &mut out);

    // 2 / 3
    let tmp = rug::Float::with_val(nd_fkind.rug_aux_prec(), 2) / 3;
    render_const(nd_fkind, "FRAC_2_3_EX", tmp, &mut out);

    // 0.4
    let tmp = rug::Float::with_val(nd_fkind.rug_aux_prec(), 4) / 10;
    render_const(nd_fkind, "FRAC_4_10_EX", tmp, &mut out);

    Ok(out)
}

pub(in super::super) fn gen_ln_special_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs) = arg_utils::parse_2_args(args)?;

    let mut out = String::new();

    let func = "(log1p(x) - log(1 - x) - 2 * x) / x";
    let poly_i = (1..=num_coeffs).map(|i| i * 2).collect::<Vec<_>>();
    let range0 = 0.1716;
    let range = (-range0, range0);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 0, "K", &mut out);

    Ok(out)
}

pub(in super::super) fn gen_ln_special_poly_ex(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs) = arg_utils::parse_2_args(args)?;

    let mut out = String::new();

    let func = "(log1p(x) - log(1 - x) - 2 * x - (2/3) * x^3 - 0.4 * x^5) / x";
    let poly_i = (1..=num_coeffs).map(|i| i * 2 + 4).collect::<Vec<_>>();
    let range = (-0.1, 0.1716);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 0, "K", &mut out);

    Ok(out)
}
