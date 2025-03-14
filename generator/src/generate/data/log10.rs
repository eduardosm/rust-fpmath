use super::super::{arg_utils, render_const, split_hi_lo, FloatKind};

pub(in super::super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let fkind: FloatKind = arg_utils::parse_1_arg(args)?;
    let aux_prec = fkind.rug_aux_prec();

    let mut out = String::new();

    // log10(e)
    let tmp = rug::Float::with_val(aux_prec, 1).exp().log10();
    let (hi, lo) = split_hi_lo(tmp, fkind.split_prec());
    render_const(fkind, "LOG10_E_HI", hi, &mut out);
    render_const(fkind, "LOG10_E_LO", lo, &mut out);

    // log10(2)
    let tmp = rug::Float::with_val(aux_prec, 2).log10();
    let (hi, lo) = split_hi_lo(tmp, fkind.split_prec());
    render_const(fkind, "LOG10_2_HI", hi, &mut out);
    render_const(fkind, "LOG10_2_LO", lo, &mut out);

    Ok(out)
}
