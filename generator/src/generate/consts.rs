use super::{FloatKind, render_const};

pub(super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let mut args = args.iter().copied();

    let Some(fkind) = args.next() else {
        return Err("not enough arguments".into());
    };
    let fkind = fkind
        .parse::<FloatKind>()
        .map_err(|_| format!("invalid aux float kind: {fkind:?}"))?;

    let rug_prec = fkind.rug_aux_prec();

    let mut out = String::new();

    for name in args {
        let value = match name {
            "INV_CBRT_2" => rug::Float::with_val(rug_prec, 2).cbrt().recip(),
            "INV_CBRT_4" => rug::Float::with_val(rug_prec, 4).cbrt().recip(),
            "LN_2" => rug::Float::with_val(rug_prec, 2).ln(),
            "LN_10" => rug::Float::with_val(rug_prec, 10).ln(),
            "LOG2_E" => rug::Float::with_val(rug_prec, 1).exp().log2(),
            "LOG2_10" => rug::Float::with_val(rug_prec, 10).log2(),
            "LOG10_E" => rug::Float::with_val(rug_prec, 1).exp().log10(),
            "LOG10_2" => rug::Float::with_val(rug_prec, 2).log10(),
            "PI" => rug::Float::with_val(rug_prec, rug::float::Constant::Pi),
            "FRAC_PI_2" => rug::Float::with_val(rug_prec, rug::float::Constant::Pi) / 2,
            "FRAC_PI_4" => rug::Float::with_val(rug_prec, rug::float::Constant::Pi) / 4,
            "FRAC_1_PI" => rug::Float::with_val(rug_prec, rug::float::Constant::Pi).recip(),
            "FRAC_3PI_4" => (rug::Float::with_val(rug_prec, rug::float::Constant::Pi) * 3) / 4,
            "FRAC_2_PI" => 2 / rug::Float::with_val(rug_prec, rug::float::Constant::Pi),
            "FRAC_PI_180" => rug::Float::with_val(rug_prec, rug::float::Constant::Pi) / 180,
            "FRAC_180_PI" => 180 / rug::Float::with_val(rug_prec, rug::float::Constant::Pi),
            "LN_PI" => rug::Float::with_val(rug_prec, rug::float::Constant::Pi).ln(),
            _ => {
                return Err(format!("unknown constant: {name:?}"));
            }
        };

        render_const(fkind, name, value, &mut out);
    }

    Ok(out)
}
