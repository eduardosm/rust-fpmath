use std::fmt::Write as _;

use super::super::arg_utils;

pub(in super::super) fn gen_frac_2_pi_large(args: &[&str]) -> Result<String, String> {
    arg_utils::expect_0_args(args)?;

    let mut out = String::new();

    let num_words = 66;
    let tmp = 2u8 / rug::Float::with_val(num_words * 24, rug::float::Constant::Pi);

    render_const_24bit_words("FRAC_2_PI_LARGE", tmp, num_words, &mut out);

    Ok(out)
}

pub(in super::super) fn gen_frac_pi_2_medium(args: &[&str]) -> Result<String, String> {
    arg_utils::expect_0_args(args)?;

    let mut out = String::new();

    let num_words = 8;
    let mut tmp = rug::Float::with_val(num_words * 24, rug::float::Constant::Pi) / 2u8;
    tmp /= 2u8;

    render_const_24bit_words("FRAC_PI_2_MEDIUM", tmp, num_words, &mut out);

    Ok(out)
}

fn render_const_24bit_words(name: &str, mut tmp: rug::Float, num_words: u32, out: &mut String) {
    out.push_str("const ");
    out.push_str(name);
    out.push_str(": &[u32] = &[\n");

    let words_per_line = 9;
    for i in 0..num_words {
        if (i % words_per_line) == 0 {
            out.push_str("    ");
        } else {
            out.push(' ');
        }
        tmp <<= 24;
        let word = tmp
            .to_integer_round(rug::float::Round::Zero)
            .unwrap()
            .0
            .to_u32()
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
}
