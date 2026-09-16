use super::super::{FloatKind, arg_utils, render_const};

pub(in super::super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let fkind: FloatKind = arg_utils::parse_1_arg(args)?;
    let sd_fkind = fkind.to_semi_double();

    let mut out = String::new();

    // π/180
    let tmp = rug::Float::with_val(sd_fkind.rug_aux_prec(), rug::float::Constant::Pi) / 180u8;
    render_const(fkind, "DEG_TO_RAD", tmp.clone(), &mut out);
    render_const(sd_fkind, "DEG_TO_RAD_EX", tmp, &mut out);

    Ok(out)
}
