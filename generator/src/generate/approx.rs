use std::fmt::Write as _;

use super::{FloatKind, arg_utils, julia, render_const_dec_value, render_const_value, sollya};

pub(super) fn gen_inv_cbrt_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs) = arg_utils::parse_2_args(args)?;

    let mut out = String::new();

    let func = "x^(-1/3)";
    let poly_i = (0..num_coeffs).collect::<Vec<_>>();
    let range = (1.0 - 0.001, 2.0 + 0.001);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 0, "K", &mut out);

    Ok(out)
}

pub(super) fn gen_exp_m1_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs, range_start, range_end) = arg_utils::parse_4_args(args)?;

    let mut out = String::new();

    let func = "expm1(x) / x - 1";
    let poly_i = (1..=num_coeffs).collect::<Vec<_>>();
    let range = (range_start, range_end);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 1, "K", &mut out);

    Ok(out)
}

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
    let (scale_fkind, fkind, bits): (FloatKind, FloatKind, u32) = arg_utils::parse_3_args(args)?;

    let mut out = String::new();

    let ftype = fkind.name();
    let prec = fkind.rug_aux_prec();
    let scale_prec = scale_fkind.float_prec();
    let num = 1 << bits;

    // The table must be consistent with `LN_LO_SCALE_TBL`, which holds the
    // rounded reciprocals, not the exact ones. Otherwise, the mismatch
    // between `ln(1 + i / num)` and `-ln(LN_LO_SCALE_TBL[i])` would introduce
    // additional error.
    writeln!(out, "// LN_TBL[i] = -ln(LN_LO_SCALE_TBL[i])").unwrap();
    writeln!(out, "static LN_TBL: [{ftype}; {}] = [", num + 1).unwrap();
    for x in 0..=num {
        let v = ((rug::Float::with_val(prec, x) >> bits) + 1u8).recip();
        // Round the reciprocal exactly as in `LN_LO_SCALE_TBL`.
        let (v, _) = rug::Float::with_val_round(scale_prec, v, rug::float::Round::Nearest);
        let ln_v = -rug::Float::with_val(prec, v).ln();
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
    writeln!(out, "static LN_LO_SCALE_TBL: [{ftype}; {}] = [", num + 1).unwrap();
    for x in 0..=num {
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

pub(super) fn gen_tan_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs) = arg_utils::parse_2_args(args)?;

    let mut out = String::new();

    let func = "tan(x) / x - 1";
    let poly_i = (1..=num_coeffs).map(|i| i * 2).collect::<Vec<_>>();
    let range0 = 0.393; // ~= π/8
    let range = (-range0, range0);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 1, "K", &mut out);

    Ok(out)
}

pub(super) fn gen_asin_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs) = arg_utils::parse_2_args(args)?;

    let mut out = String::new();

    let func = "asin(x) - x";
    let poly_i = (1..=num_coeffs).map(|i| i * 2 + 1).collect::<Vec<_>>();
    let range = (-0.00001, 0.50001);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 0, "K", &mut out);

    Ok(out)
}

pub(super) fn gen_atan_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs, range_start, range_end) = arg_utils::parse_4_args(args)?;

    let mut out = String::new();

    let func = "atan(x) - x";
    let poly_i = (1..=num_coeffs).map(|i| i * 2 + 1).collect::<Vec<_>>();
    let range = (range_start, range_end);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 0, "K", &mut out);

    Ok(out)
}

pub(super) fn gen_asinh_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs, range_start, range_end) = arg_utils::parse_4_args(args)?;

    let mut out = String::new();

    let func = "asinh(x) / x - 1";
    let poly_i = (1..=num_coeffs).map(|i| i * 2).collect::<Vec<_>>();
    let range = (range_start, range_end);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 1, "K", &mut out);

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
