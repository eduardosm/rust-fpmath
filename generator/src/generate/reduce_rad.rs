use std::fmt::Write as _;

use super::arg_utils;

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
