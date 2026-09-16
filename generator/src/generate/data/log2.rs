use super::super::{FloatKind, arg_utils, render_double_const, split_hi_lo};

pub(in super::super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let fkind: FloatKind = arg_utils::parse_1_arg(args)?;
    let aux_prec = fkind.rug_aux_prec();

    let mut out = String::new();

    // log2(e)
    let tmp = rug::Float::with_val(aux_prec, 1).exp().log2();
    let (hi, lo) = split_hi_lo(tmp, fkind.split_prec());
    render_double_const(fkind, "SemiDouble", "LOG2_E_EX", hi, lo, &mut out);

    Ok(out)
}
