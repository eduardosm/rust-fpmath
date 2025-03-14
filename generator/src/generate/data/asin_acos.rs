use super::super::{arg_utils, render_const, sollya, split_hi_lo, FloatKind};

pub(in super::super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let fkind: FloatKind = arg_utils::parse_1_arg(args)?;
    let aux_prec = fkind.rug_aux_prec();

    let mut out = String::new();

    // π/2
    let tmp = rug::Float::with_val(aux_prec, rug::float::Constant::Pi) / 2u8;
    let (hi, lo) = split_hi_lo(tmp, fkind.float_prec());
    render_const(fkind, "FRAC_PI_2_HI", hi, &mut out);
    render_const(fkind, "FRAC_PI_2_LO", lo, &mut out);

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
