use std::fmt::Write as _;

use super::{FloatKind, arg_utils, render_const};

pub(super) fn gen_medium_consts(args: &[&str]) -> Result<String, String> {
    let (fkind, fkind_ex): (FloatKind, FloatKind) = arg_utils::parse_2_args(args)?;

    let mut out = String::new();

    // π/2
    let tmp = rug::Float::with_val(
        fkind.rug_aux_prec() + fkind_ex.rug_aux_prec(),
        rug::float::Constant::Pi,
    ) / 2;

    let (p0base, p0ex) = split_ex(tmp, fkind.float_prec());
    render_const(fkind, "FRAC_PI_2_P0", p0base, &mut out);
    render_const(fkind_ex, "FRAC_PI_2_P0EX", p0ex.clone(), &mut out);

    Ok(out)
}

fn split_ex(tmp: rug::Float, base_prec: u32) -> (rug::Float, rug::Float) {
    let (base, _) = rug::Float::with_val_round(base_prec, &tmp, rug::float::Round::Zero);
    let ex = tmp - &base;
    (base, ex)
}

pub(super) fn gen_frac_2_pi_large(args: &[&str]) -> Result<String, String> {
    arg_utils::expect_0_args(args)?;

    let mut out = String::new();

    let num_words = 66;
    let tmp = 2u8 / rug::Float::with_val(num_words * 24, rug::float::Constant::Pi);

    render_const_24bit_words("FRAC_2_PI_LARGE", tmp, num_words, &mut out);

    Ok(out)
}

pub(super) fn gen_frac_pi_2_medium(args: &[&str]) -> Result<String, String> {
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
    if num_words % words_per_line != 0 {
        out.push('\n');
    }
    out.push_str("];\n");
}
