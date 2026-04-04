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
            "PI" => rug::Float::with_val(rug_prec, rug::float::Constant::Pi),
            "FRAC_PI_2" => rug::Float::with_val(rug_prec, rug::float::Constant::Pi) / 2,
            "FRAC_PI_4" => rug::Float::with_val(rug_prec, rug::float::Constant::Pi) / 4,
            "FRAC_2_PI" => 2 / rug::Float::with_val(rug_prec, rug::float::Constant::Pi),
            _ => {
                return Err(format!("unknown constant: {name:?}"));
            }
        };

        render_const(fkind, name, value, &mut out);
    }

    Ok(out)
}
