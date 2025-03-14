use super::super::{arg_utils, render_const, sollya, split_hi_lo, FloatKind};

pub(in super::super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let fkind: FloatKind = arg_utils::parse_1_arg(args)?;
    let aux_prec = fkind.rug_aux_prec();

    let mut out = String::new();

    // cbrt(2)
    let tmp = rug::Float::with_val(aux_prec, 2).cbrt();
    let (hi, lo) = split_hi_lo(tmp, fkind.split_prec());
    render_const(fkind, "CBRT_2_HI", hi, &mut out);
    render_const(fkind, "CBRT_2_LO", lo, &mut out);

    // cbrt(4)
    let tmp = rug::Float::with_val(aux_prec, 4).cbrt();
    let (hi, lo) = split_hi_lo(tmp, fkind.split_prec());
    render_const(fkind, "CBRT_4_HI", hi, &mut out);
    render_const(fkind, "CBRT_4_LO", lo, &mut out);

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
