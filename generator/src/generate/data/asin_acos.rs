use super::super::{FloatKind, arg_utils, render_const, sollya};

pub(in super::super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let fkind: FloatKind = arg_utils::parse_1_arg(args)?;
    let nd_fkind = fkind.to_double();

    let mut out = String::new();

    // π/2
    let tmp = rug::Float::with_val(nd_fkind.rug_aux_prec(), rug::float::Constant::Pi) / 2u8;
    render_const(nd_fkind, "FRAC_PI_2_EX", tmp, &mut out);

    Ok(out)
}

pub(in super::super) fn gen_asin_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs) = arg_utils::parse_2_args(args)?;

    let mut out = String::new();

    let func = "asin(x) - x";
    let poly_i = (1..=num_coeffs).map(|i| i * 2 + 1).collect::<Vec<_>>();
    let range = (-0.001, 0.501);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, -3, "K", &mut out);

    Ok(out)
}
