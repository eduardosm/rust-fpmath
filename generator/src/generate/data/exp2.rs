use super::super::{arg_utils, render_const, FloatKind};

pub(in super::super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let fkind: FloatKind = arg_utils::parse_1_arg(args)?;
    let aux_prec = fkind.rug_aux_prec();

    let mut out = String::new();

    // ln(2)
    let tmp = rug::Float::with_val(aux_prec, rug::float::Constant::Log2);
    render_const(fkind, "LN_2", tmp, &mut out);

    Ok(out)
}
