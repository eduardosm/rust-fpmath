use super::super::{FloatKind, arg_utils, render_const};

pub(in super::super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let fkind: FloatKind = arg_utils::parse_1_arg(args)?;
    let sd_fkind = fkind.to_semi_double();

    let mut out = String::new();

    // 1/π
    let tmp = 1u8 / rug::Float::with_val(sd_fkind.rug_aux_prec(), rug::float::Constant::Pi);
    render_const(sd_fkind, "FRAC_1_PI_EX", tmp, &mut out);

    Ok(out)
}
