use super::super::{arg_utils, render_const, split_hi_lo, FloatKind};

pub(in super::super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let fkind: FloatKind = arg_utils::parse_1_arg(args)?;
    let aux_prec = fkind.rug_aux_prec();

    let mut out = String::new();

    // 1/π
    let tmp = 1u8 / rug::Float::with_val(aux_prec, rug::float::Constant::Pi);
    let (hi, lo) = split_hi_lo(tmp, fkind.split_prec());
    render_const(fkind, "FRAC_1_PI_HI", hi, &mut out);
    render_const(fkind, "FRAC_1_PI_LO", lo, &mut out);

    Ok(out)
}
