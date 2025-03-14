use super::super::{arg_utils, render_const, split_hi_lo, FloatKind};

pub(in super::super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let fkind: FloatKind = arg_utils::parse_1_arg(args)?;
    let aux_prec = fkind.rug_aux_prec();

    let mut out = String::new();

    // 180/π
    let tmp = 180u8 / rug::Float::with_val(aux_prec, rug::float::Constant::Pi);
    let (hi, lo) = split_hi_lo(tmp.clone(), fkind.split_prec());
    render_const(fkind, "RAD_TO_DEG", tmp, &mut out);
    render_const(fkind, "RAD_TO_DEG_HI", hi, &mut out);
    render_const(fkind, "RAD_TO_DEG_LO", lo, &mut out);

    Ok(out)
}
