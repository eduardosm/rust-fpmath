use std::fmt::Write as _;

use super::{FloatKind, arg_utils, render_const};

/// Splits `π * 2^e` into `f64` parts `{NAME}0`, `{NAME}1`, ..., where each
/// part `i` has at most `bits[i]` significant bits. All parts except the
/// last one are truncated (rounded towards zero), and the last one is
/// rounded to nearest.
///
/// Arguments: `NAME e bits...`
pub(super) fn gen_split_pi(args: &[&str]) -> Result<String, String> {
    let [name, e, bits @ ..] = args else {
        return Err(format!(
            "expected at least 3 arguments, found {}",
            args.len()
        ));
    };
    if bits.is_empty() {
        return Err("expected at least 3 arguments, found 2".into());
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
    let mut rem = rug::Float::with_val(prec, rug::float::Constant::Pi);
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

/// Bits of `2^-shift / π`, in 64-bit words, most significant first:
/// `FRAC_1_PI_BITS[i]` holds bits `64*i + 1 ..= 64*i + 64` of the binary
/// expansion (i.e., `2^-shift / π ~= sum(FRAC_1_PI_BITS[i] * 2^(-64 * (i + 1)))`).
///
/// Arguments: `shift num_words`
pub(super) fn gen_frac_1_pi_bits(args: &[&str]) -> Result<String, String> {
    let (shift, num_words): (u32, u32) = arg_utils::parse_2_args(args)?;

    let mut out = String::new();

    let prec = (num_words + 2) * 64 + shift;
    let mut tmp = rug::Float::with_val(prec, rug::float::Constant::Pi).recip();
    tmp >>= shift;

    writeln!(
        out,
        "// 2^-{shift} / π ~= sum(FRAC_1_PI_BITS[i] * 2^(-64 * (i + 1)))"
    )
    .unwrap();
    writeln!(out, "static FRAC_1_PI_BITS: [u64; {num_words}] = [").unwrap();
    for _ in 0..num_words {
        tmp <<= 64;
        let word = tmp
            .to_integer_round(rug::float::Round::Zero)
            .unwrap()
            .0
            .to_u64()
            .unwrap();
        writeln!(out, "    0x{word:016X},").unwrap();
        tmp -= word;
    }
    out.push_str("];\n");

    Ok(out)
}
