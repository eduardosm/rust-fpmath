use std::fmt::Write as _;

use super::super::arg_utils;

pub(in super::super) fn gen_frac_2_pi_large(args: &[&str]) -> Result<String, String> {
    arg_utils::expect_0_args(args)?;

    let mut out = String::new();

    let num_words = 66;
    let mut tmp = 2u8 / rug::Float::with_val(num_words * 24, rug::float::Constant::Pi);

    out.push_str("const FRAC_2_PI_LARGE: &[u32] = &[\n");
    let words_per_line = 9;
    for i in 0..num_words {
        if (i % words_per_line) == 0 {
            out.push_str("    ");
        } else {
            out.push(' ');
        }
        tmp <<= 24;
        let word = tmp
            .to_u32_saturating_round(rug::float::Round::Zero)
            .unwrap();
        write!(out, "0x{word:06X},").unwrap();
        tmp -= word;
        if (i % words_per_line) == (words_per_line - 1) {
            out.push('\n');
        }
    }
    if (num_words % words_per_line) != 0 {
        out.push('\n');
    }
    out.push_str("];\n");

    Ok(out)
}

pub(in super::super) fn gen_frac_pi_2_medium(args: &[&str]) -> Result<String, String> {
    arg_utils::expect_0_args(args)?;

    let mut out = String::new();

    let num_words = 8;
    let mut tmp = rug::Float::with_val(num_words * 24, rug::float::Constant::Pi) / 2u8;

    out.push_str("const FRAC_PI_2_MEDIUM: &[u32] = &[\n");
    tmp /= 2u8;
    let words_per_line = 9;
    for i in 0..num_words {
        if (i % words_per_line) == 0 {
            out.push_str("    ");
        } else {
            out.push(' ');
        }
        tmp <<= 24;
        let word = tmp
            .to_u32_saturating_round(rug::float::Round::Zero)
            .unwrap();
        write!(out, "0x{word:06X},").unwrap();
        tmp -= word;
        if (i % words_per_line) == (words_per_line - 1) {
            out.push('\n');
        }
    }
    if (num_words % words_per_line) != 0 {
        out.push('\n');
    }
    out.push_str("];\n");

    Ok(out)
}
