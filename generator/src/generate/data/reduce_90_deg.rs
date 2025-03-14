use super::super::{arg_utils, render_const, split_hi_lo, FloatKind};

pub(in super::super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let fkind: FloatKind = arg_utils::parse_1_arg(args)?;
    let aux_prec = fkind.rug_aux_prec();

    let mut out = String::new();

    // π/180
    let tmp = rug::Float::with_val(aux_prec, rug::float::Constant::Pi) / 180u8;
    let (hi, lo) = split_hi_lo(tmp.clone(), fkind.split_prec());
    render_const(fkind, "DEG_TO_RAD", tmp, &mut out);
    render_const(fkind, "DEG_TO_RAD_HI", hi, &mut out);
    render_const(fkind, "DEG_TO_RAD_LO", lo, &mut out);

    Ok(out)
}
