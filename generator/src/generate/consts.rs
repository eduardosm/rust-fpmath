use super::{FloatKind, render_const};

/// Generates the named mathematical constants (e.g., `PI`, `LN_2`) with type
/// `fkind`.
///
/// Arguments: `fkind NAME...`
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
        let value = named_const(name, rug_prec)?;
        render_const(fkind, name, value, &mut out);
    }

    Ok(out)
}

/// Splits the named mathematical constant `CONST` (see `consts`) multiplied
/// by `2^e` into `f64` parts `{NAME}0`, `{NAME}1`, ..., where each part `i`
/// has at most `bits[i]` significant bits. All parts except the last one are
/// truncated (rounded towards zero), and the last one is rounded to nearest.
///
/// Arguments: `NAME CONST e bits...`
pub(super) fn gen_split_const(args: &[&str]) -> Result<String, String> {
    let [name, const_name, e, bits @ ..] = args else {
        return Err(format!(
            "expected at least 4 arguments, found {}",
            args.len()
        ));
    };
    if bits.is_empty() {
        return Err("expected at least 4 arguments, found 3".into());
    }
    let e: i32 = e
        .parse()
        .map_err(|e| format!("failed to parse exponent: {e}"))?;
    let bits = bits
        .iter()
        .map(|b| match b.parse::<u32>() {
            Ok(b @ 1..=53) => Ok(b),
            Ok(b) => Err(format!("invalid number of bits: {b}")),
            Err(e) => Err(format!("failed to parse number of bits {b:?}: {e}")),
        })
        .collect::<Result<Vec<u32>, _>>()?;

    let mut out = String::new();

    let prec = 1024;
    let mut rem = named_const(const_name, prec)?;
    if e >= 0 {
        rem <<= e as u32;
    } else {
        rem >>= e.unsigned_abs();
    }

    for (i, &part_bits) in bits.iter().enumerate() {
        let round = if i == bits.len() - 1 {
            rug::float::Round::Nearest
        } else {
            rug::float::Round::Zero
        };
        let (part, _) = rug::Float::with_val_round(part_bits, &rem, round);
        rem -= &part;
        render_const(FloatKind::F64, &format!("{name}{i}"), part, &mut out);
    }

    Ok(out)
}

/// Returns the value of the named mathematical constant `name` (e.g., `PI`,
/// `LN_2`) with precision `prec`.
fn named_const(name: &str, prec: u32) -> Result<rug::Float, String> {
    let value = match name {
        "CBRT_2" => rug::Float::with_val(prec, 2).cbrt(),
        "CBRT_4" => rug::Float::with_val(prec, 4).cbrt(),
        "LN_2" => rug::Float::with_val(prec, 2).ln(),
        "LN_10" => rug::Float::with_val(prec, 10).ln(),
        "LOG2_E" => rug::Float::with_val(prec, 1).exp().log2(),
        "LOG2_10" => rug::Float::with_val(prec, 10).log2(),
        "LOG10_E" => rug::Float::with_val(prec, 1).exp().log10(),
        "LOG10_2" => rug::Float::with_val(prec, 2).log10(),
        "PI" => rug::Float::with_val(prec, rug::float::Constant::Pi),
        "FRAC_PI_2" => rug::Float::with_val(prec, rug::float::Constant::Pi) / 2,
        "FRAC_PI_4" => rug::Float::with_val(prec, rug::float::Constant::Pi) / 4,
        "FRAC_1_PI" => rug::Float::with_val(prec, rug::float::Constant::Pi).recip(),
        "FRAC_3PI_4" => (rug::Float::with_val(prec, rug::float::Constant::Pi) * 3) / 4,
        "FRAC_2_PI" => 2 / rug::Float::with_val(prec, rug::float::Constant::Pi),
        "FRAC_PI_180" => rug::Float::with_val(prec, rug::float::Constant::Pi) / 180,
        "FRAC_180_PI" => 180 / rug::Float::with_val(prec, rug::float::Constant::Pi),
        "LN_PI" => rug::Float::with_val(prec, rug::float::Constant::Pi).ln(),
        _ => {
            return Err(format!("unknown constant: {name:?}"));
        }
    };
    Ok(value)
}
