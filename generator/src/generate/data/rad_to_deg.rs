use super::super::{FloatKind, arg_utils, render_const};

pub(in super::super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let fkind: FloatKind = arg_utils::parse_1_arg(args)?;
    let sd_fkind = fkind.to_semi_double();

    let mut out = String::new();

    // 180/π
    let tmp = 180u8 / rug::Float::with_val(sd_fkind.rug_aux_prec(), rug::float::Constant::Pi);
    render_const(fkind, "RAD_TO_DEG", tmp.clone(), &mut out);
    render_const(sd_fkind, "RAD_TO_DEG_EX", tmp, &mut out);

    Ok(out)
}
