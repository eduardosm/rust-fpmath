use std::fmt::Write as _;

use super::{FloatKind, arg_utils, julia, render_const_dec_value, render_const_value, sollya};

pub(super) fn gen_ln_1p_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs, range_start, range_end) = arg_utils::parse_4_args(args)?;

    let mut out = String::new();

    let func = "log1p(x) / x - 1";
    let poly_i = (1..=num_coeffs).collect::<Vec<_>>();
    let range = (range_start, range_end);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 1, "K", &mut out);

    Ok(out)
}

pub(super) fn gen_ln_table(args: &[&str]) -> Result<String, String> {
    let (fkind, bits): (FloatKind, u32) = arg_utils::parse_2_args(args)?;

    let mut out = String::new();

    let ftype = fkind.name();
    let prec = fkind.rug_aux_prec();
    let num = 1 << bits;

    writeln!(
        out,
        "// LN_TBL[i] = ln(1 + i / {num})       if i < {}",
        num / 2,
    )
    .unwrap();
    writeln!(
        out,
        "//           = ln((1 + i / {num}) / 2) if i >= {}",
        num / 2,
    )
    .unwrap();
    writeln!(out, "static LN_TBL: [{ftype}; {num}] = [").unwrap();
    for x in 0..num {
        let mut v = (rug::Float::with_val(prec, x) >> bits) + 1u8;
        if x >= num / 2 {
            // The upper half of the table contains `ln(v / 2)` instead of `ln(v)`
            // to avoid cancellation when `x` is close to 1.
            v >>= 1;
        }
        let ln_v = v.ln();
        out.push_str("    ");
        render_const_value(fkind, &ln_v, &mut out);
        out.push_str(", // ");
        render_const_dec_value(fkind, &ln_v, &mut out);
        out.push('\n');
    }
    writeln!(out, "];").unwrap();

    Ok(out)
}

pub(super) fn gen_ln_lo_scale_table(args: &[&str]) -> Result<String, String> {
    let (fkind, bits): (FloatKind, u32) = arg_utils::parse_2_args(args)?;

    let mut out = String::new();

    let ftype = fkind.name();
    let prec = fkind.rug_aux_prec();
    let num = 1 << bits;

    writeln!(out, "// LN_LO_SCALE_TBL[i] = 1 / (1 + i / {num})").unwrap();
    writeln!(out, "static LN_LO_SCALE_TBL: [{ftype}; {num}] = [").unwrap();
    for x in 0..num {
        let v = (rug::Float::with_val(prec, x) >> bits) + 1u8;
        let v = v.recip();
        out.push_str("    ");
        render_const_value(fkind, &v, &mut out);
        out.push_str(", // ");
        render_const_dec_value(fkind, &v, &mut out);
        out.push('\n');
    }
    writeln!(out, "];").unwrap();

    Ok(out)
}

pub(super) fn gen_sin_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs) = arg_utils::parse_2_args(args)?;

    let mut out = String::new();

    let func = "sin(x) / x - 1";
    let poly_i = (1..=num_coeffs).map(|i| i * 2).collect::<Vec<_>>();
    let range0 = 0.786; // ~= π/4
    let range = (-range0, range0);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 1, "K", &mut out);

    Ok(out)
}

pub(super) fn gen_cos_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs) = arg_utils::parse_2_args(args)?;

    let mut out = String::new();

    let func = "cos(x) - (1 - 0.5 * x^2)";
    let poly_i = (1..=num_coeffs).map(|i| i * 2 + 2).collect::<Vec<_>>();
    let range0 = 0.786; // ~= π/4
    let range = (-range0, range0);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 0, "K", &mut out);

    Ok(out)
}

pub(super) fn gen_gamma_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, poly_deg, offset, range_start, range_end): (_, i32, f64, f64, f64) =
        arg_utils::parse_5_args(args)?;

    let mut out = String::new();

    let func = format!("SpecialFunctions.gamma(x + {offset})");
    let wfunc = "1 / fx";
    let range = (range_start, range_end);

    julia::run_and_render_remez(fkind, &func, wfunc, range, poly_deg, 0, "K", &mut out);

    Ok(out)
}

pub(super) fn gen_ln_gamma_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, poly_deg, offset, range_start, range_end): (_, i32, f64, f64, f64) =
        arg_utils::parse_5_args(args)?;

    let mut out = String::new();

    let func = format!("SpecialFunctions.lgamma(x + {offset}) / x");
    let wfunc = "1 / fx";
    let range = (range_start, range_end);

    julia::run_and_render_remez(fkind, &func, wfunc, range, poly_deg - 1, 1, "K", &mut out);

    Ok(out)
}

pub(super) fn gen_gamma_lanczos_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, poly_deg, g, range_start, range_end): (_, i32, f64, f64, f64) =
        arg_utils::parse_5_args(args)?;

    let mut out = String::new();

    let g = format!("BigFloat({g})");

    // SpecialFunctions.gamma(1/x) / ((1/x + g - 0.5)^(1/x - 0.5) * exp(-(1/x + g - 0.5)))
    let func = format!(
        "exp(SpecialFunctions.lgamma(1/x) - (1/x - 0.5) * log(1/x + {g} - 0.5) + (1/x + {g} - 0.5))"
    );
    let wfunc = "1";
    let range = (range_start, range_end);

    julia::run_and_render_remez(fkind, &func, wfunc, range, poly_deg, 0, "K", &mut out);

    Ok(out)
}
