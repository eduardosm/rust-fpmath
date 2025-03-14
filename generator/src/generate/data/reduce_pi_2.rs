use super::super::{arg_utils, render_const, split_hi_lo, FloatKind};

pub(in super::super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let fkind: FloatKind = arg_utils::parse_1_arg(args)?;
    let aux_prec = fkind.rug_aux_prec();

    let mut out = String::new();

    let split_prec = match fkind {
        FloatKind::F32 => 15,
        FloatKind::F64 => 33,
    };

    // π/2
    let tmp = rug::Float::with_val(aux_prec * 3, rug::float::Constant::Pi) / 2;
    let (hi, hiex) = split_hi_lo(tmp, split_prec);
    let (mi, miex) = split_hi_lo(hiex.clone(), split_prec);
    let (lo, loex) = split_hi_lo(miex.clone(), split_prec);
    render_const(fkind, "FRAC_PI_2_HI", hi, &mut out);
    render_const(fkind, "FRAC_PI_2_HIEX", hiex, &mut out);
    render_const(fkind, "FRAC_PI_2_MI", mi, &mut out);
    render_const(fkind, "FRAC_PI_2_MIEX", miex, &mut out);
    render_const(fkind, "FRAC_PI_2_LO", lo, &mut out);
    render_const(fkind, "FRAC_PI_2_LOEX", loex, &mut out);

    Ok(out)
}
