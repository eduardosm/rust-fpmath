use super::super::{FloatKind, arg_utils, render_const};

pub(in super::super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let fkind: FloatKind = arg_utils::parse_1_arg(args)?;
    let sd_fkind = fkind.to_semi_double();

    let mut out = String::new();

    // log2(e)
    let tmp = rug::Float::with_val(sd_fkind.rug_aux_prec(), 1)
        .exp()
        .log2();
    render_const(sd_fkind, "LOG2_E_EX", tmp, &mut out);

    Ok(out)
}
