use super::super::{FloatKind, arg_utils, render_const, sollya};

pub(in super::super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let fkind: FloatKind = arg_utils::parse_1_arg(args)?;
    let sd_fkind = fkind.to_semi_double();

    let mut out = String::new();

    // cbrt(2)
    let tmp = rug::Float::with_val(sd_fkind.rug_aux_prec(), 2).cbrt();
    render_const(sd_fkind, "CBRT_2_EX", tmp, &mut out);

    // cbrt(4)
    let tmp = rug::Float::with_val(sd_fkind.rug_aux_prec(), 4).cbrt();
    render_const(sd_fkind, "CBRT_4_EX", tmp, &mut out);

    Ok(out)
}

pub(in super::super) fn gen_inv_cbrt_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs) = arg_utils::parse_2_args(args)?;

    let mut out = String::new();

    let func = "x^(-1/3)";
    let poly_i = (0..num_coeffs).collect::<Vec<_>>();
    let range = (1.0 - 0.001, 2.0 + 0.001);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 0, "K", &mut out);

    Ok(out)
}
