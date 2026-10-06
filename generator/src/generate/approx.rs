use std::fmt::Write as _;

use super::{FloatKind, arg_utils, julia, render_const_dec_value, render_const_value, sollya};

/// Generates a table of 128 approximations of `1 / sqrt(m)` for `m` in
/// `[1, 4)`, in 0.16 fixed point. Entry `i` covers the subinterval
/// `[a, b]`, with `a = 2^(i >> 6) * (1 + (i & 63) / 64)` and
/// `b = 2^(i >> 6) * (1 + ((i & 63) + 1) / 64)`, and minimizes the maximum
/// relative error in it.
///
/// Arguments: none
pub(super) fn gen_rsqrt_table(args: &[&str]) -> Result<String, String> {
    arg_utils::expect_0_args(args)?;

    let mut out = String::new();

    let prec = 128;

    writeln!(
        out,
        "// RSQRT_TBL[i] ~= 1 / sqrt(m) (0.16 fixed point), minimizing the maximum"
    )
    .unwrap();
    writeln!(out, "// relative error for m in [a, b], with").unwrap();
    writeln!(out, "// * a = 2^(i >> 6) * (1 + (i & 63) / 64)").unwrap();
    writeln!(out, "// * b = 2^(i >> 6) * (1 + ((i & 63) + 1) / 64)").unwrap();
    writeln!(out, "static RSQRT_TBL: [u16; 128] = [").unwrap();
    for i in 0..128u32 {
        let a = rug::Float::with_val(prec, 64 + (i & 63)) << (i >> 6);
        let b = rug::Float::with_val(prec, 65 + (i & 63)) << (i >> 6);
        // r = 2 / (sqrt(a) + sqrt(b)) minimizes max(|r * sqrt(m) - 1|) for m in [a, b]
        // (a and b are scaled by 64, so sqrt(64) = 8 is compensated)
        let r = rug::Float::with_val(prec, a.sqrt() + b.sqrt()).recip() * 16u8;
        let v = (r << 16u32).to_integer().unwrap().to_u16().unwrap();
        writeln!(out, "    0x{v:04X},").unwrap();
    }
    writeln!(out, "];").unwrap();

    Ok(out)
}

/// Generates an approximation of `x^(1/3)` in `[1 - 0.001, 2 + 0.001]` with
/// the powers `1, x, x^2, ...`.
///
/// Arguments: `fkind num_coeffs`
pub(super) fn gen_cbrt_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs) = arg_utils::parse_2_args(args)?;

    let mut out = String::new();

    let func = "x^(1/3)";
    let poly_i = (0..num_coeffs).collect::<Vec<_>>();
    let range = (1.0 - 0.001, 2.0 + 0.001);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 0, "K", &mut out);

    Ok(out)
}

/// Generates an approximation of `((1 + x)^(1/3) - 1) / x` in
/// `[range_start, range_end]` with the powers `1, x, x^2, ...`.
///
/// Arguments: `fkind num_coeffs range_start range_end`
pub(super) fn gen_cbrt_1p_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs, range_start, range_end) = arg_utils::parse_4_args(args)?;

    let mut out = String::new();

    let func = "((1 + x)^(1/3) - 1) / x";
    let poly_i = (0..num_coeffs).collect::<Vec<_>>();
    let range = (range_start, range_end);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 0, "K", &mut out);

    Ok(out)
}

/// Generates an approximation of `(exp(x) - 1) / x - 1` in
/// `[range_start, range_end]` with the powers `x, x^2, ...`.
///
/// When `num_fixed` (default 1) is greater than 1, the coefficients of
/// `x^1, ..., x^(num_fixed - 1)` are fixed to those of the Taylor series
/// (`1/2!, 1/3!, ...`), so it approximates
/// `(exp(x) - 1) / x - 1 - x / 2! - ... - x^(num_fixed - 1) / num_fixed!` with
/// the powers `x^num_fixed, x^(num_fixed + 1), ...`.
///
/// Arguments: `fkind num_coeffs range_start range_end [num_fixed]`
pub(super) fn gen_exp_m1_poly(args: &[&str]) -> Result<String, String> {
    let (args, num_fixed) = match args {
        [args @ .., num_fixed] if args.len() == 4 => {
            let num_fixed: i32 = num_fixed
                .parse()
                .map_err(|e| format!("failed to parse fifth argument {num_fixed:?}: {e}"))?;
            if num_fixed < 1 {
                return Err(format!("invalid number of fixed terms: {num_fixed}"));
            }
            (args, num_fixed)
        }
        [_, _, _, _] => (args, 1),
        _ => return Err(format!("expected 4 or 5 arguments, found {}", args.len())),
    };
    let (fkind, num_coeffs, range_start, range_end): (_, i32, _, _) =
        arg_utils::parse_4_args(args)?;

    let mut out = String::new();

    let mut func = String::from("expm1(x) / x - 1");
    let mut factorial = 1u64;
    for i in 2..=num_fixed {
        factorial = factorial.checked_mul(i as u64).unwrap();
        write!(func, " - x^{} / {factorial}", i - 1).unwrap();
    }
    let poly_i = (num_fixed..(num_fixed + num_coeffs)).collect::<Vec<_>>();
    let range = (range_start, range_end);

    sollya::run_and_render_remez(fkind, &func, range, &poly_i, 1, "K", &mut out);

    Ok(out)
}

/// Generates a table of `2^(i / 2^bits)` for `i` in `0..2^bits`. Each entry
/// is a pair `(T, υ)` (as the bits of `f64` values), where `T` is
/// `2^(i / 2^bits)` rounded to `f64` and `υ` is `(2^(i / 2^bits) - Th) / T`
/// rounded to `f64`, with `Th` being `T` truncated to `hi_bits` significant
/// bits (so `2^(i / 2^bits) ~= Th + T * υ`).
///
/// Arguments: `bits hi_bits`
pub(super) fn gen_exp2_table(args: &[&str]) -> Result<String, String> {
    let (bits, hi_bits): (u32, u32) = arg_utils::parse_2_args(args)?;
    if !(1..=10).contains(&bits) {
        return Err(format!("invalid number of bits: {bits}"));
    }
    if !(1..=53).contains(&hi_bits) {
        return Err(format!("invalid number of high bits: {hi_bits}"));
    }

    let mut out = String::new();

    let fkind = FloatKind::F64;
    let prec = fkind.rug_aux_prec();
    let num = 1u32 << bits;

    writeln!(
        out,
        "// EXP2_TBL[i] = (bits(T), bits(υ)), where T is 2^(i / {num}) rounded,"
    )
    .unwrap();
    writeln!(
        out,
        "// Th is T truncated to {hi_bits} bits and 2^(i / {num}) = Th + T * υ"
    )
    .unwrap();
    writeln!(out, "static EXP2_TBL: [(u64, u64); {num}] = [").unwrap();
    for i in 0..num {
        let v = (rug::Float::with_val(prec, i) >> bits).exp2();
        let t = rug::Float::with_val(fkind.float_prec(), &v);
        let (th, _) = rug::Float::with_val_round(hi_bits, &t, rug::float::Round::Zero);
        let upsilon = rug::Float::with_val(prec, &v - &th) / &t;
        write!(
            out,
            "    (0x{:016X}, 0x{:016X}), // ",
            t.to_f64().to_bits(),
            upsilon.to_f64().to_bits(),
        )
        .unwrap();
        render_const_dec_value(fkind, &t, &mut out);
        out.push('\n');
    }
    writeln!(out, "];").unwrap();

    Ok(out)
}

/// Generates an approximation of `ln(1 + x) / x - 1` in
/// `[range_start, range_end]` with the powers `x, x^2, ...`.
///
/// Arguments: `fkind num_coeffs range_start range_end`
pub(super) fn gen_ln_1p_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs, range_start, range_end) = arg_utils::parse_4_args(args)?;

    let mut out = String::new();

    let func = "log1p(x) / x - 1";
    let poly_i = (1..=num_coeffs).collect::<Vec<_>>();
    let range = (range_start, range_end);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 1, "K", &mut out);

    Ok(out)
}

/// Generates a table of `-ln(LN_LO_SCALE_TBL[i])` (i.e., approximately
/// `ln(1 + i / 2^bits)`) for `i` in `0..=2^bits`, with values of type
/// `fkind`.
///
/// `LN_LO_SCALE_TBL[i]` is rounded to `scale_fkind`, so `scale_fkind` and
/// `bits` must match the arguments of `ln_lo_scale_table`.
///
/// Arguments: `scale_fkind fkind bits`
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

/// Generates a table of `1 / (1 + i / 2^bits)` for `i` in `0..=2^bits`.
///
/// Arguments: `fkind bits`
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

/// Generates an approximation of `sin(x) / x - 1` in `[-range, range]` with
/// the even powers `x^2, x^4, ...`.
///
/// Arguments: `fkind num_coeffs range`
pub(super) fn gen_sin_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs, range0): (_, i32, f64) = arg_utils::parse_3_args(args)?;

    let mut out = String::new();

    let func = "sin(x) / x - 1";
    let poly_i = (1..=num_coeffs).map(|i| i * 2).collect::<Vec<_>>();
    let range = (-range0, range0);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 1, "K", &mut out);

    Ok(out)
}

/// Generates an approximation of `cos(x) - (1 - x^2 / 2)` in `[-range, range]`
/// with the even powers `x^4, x^6, ...`.
///
/// Arguments: `fkind num_coeffs range`
pub(super) fn gen_cos_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs, range0): (_, i32, f64) = arg_utils::parse_3_args(args)?;

    let mut out = String::new();

    let func = "cos(x) - (1 - 0.5 * x^2)";
    let poly_i = (1..=num_coeffs).map(|i| i * 2 + 2).collect::<Vec<_>>();
    let range = (-range0, range0);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 0, "K", &mut out);

    Ok(out)
}

/// Generates a table of `sin(i * π / n)` for `i` in `0..(2 * n)` (a full period).
///
/// The values are computed from the first quadrant, so the table is exactly
/// symmetric: `SIN_PI_TBL[2n - i] == -SIN_PI_TBL[i]` and
/// `SIN_PI_TBL[n/2 - i] == SIN_PI_TBL[n/2 + i]`.
pub(super) fn gen_sin_pi_table(args: &[&str]) -> Result<String, String> {
    let (fkind, n): (FloatKind, u32) = arg_utils::parse_2_args(args)?;
    if n == 0 || n % 2 != 0 {
        return Err(format!("n must be even and non-zero, found {n}"));
    }

    let mut out = String::new();

    let ftype = fkind.name();
    let prec = fkind.rug_aux_prec();
    let len = 2 * n;

    writeln!(out, "// SIN_PI_TBL[i] = sin(i * π / {n})").unwrap();
    writeln!(out, "static SIN_PI_TBL: [{ftype}; {len}] = [").unwrap();
    for i in 0..len {
        // sin(i * π / n) = sign * sin(k * π / n), with 0 <= k <= n / 2
        let negative = i >= n;
        let mut k = i % n;
        if k > n / 2 {
            k = n - k;
        }

        let v = if k == 0 {
            rug::Float::with_val(prec, 0)
        } else if k == n / 2 {
            rug::Float::with_val(prec, 1)
        } else {
            let angle = rug::Float::with_val(prec, rug::float::Constant::Pi) * k / n;
            angle.sin()
        };
        let v = if negative && k != 0 { -v } else { v };

        out.push_str("    ");
        render_const_value(fkind, &v, &mut out);
        out.push_str(", // ");
        render_const_dec_value(fkind, &v, &mut out);
        out.push('\n');
    }
    writeln!(out, "];").unwrap();

    Ok(out)
}

/// Generates an approximation of `tan(x) / x - 1` in `[-range, range]` with
/// the even powers `x^2, x^4, ...`.
///
/// Arguments: `fkind num_coeffs range`
pub(super) fn gen_tan_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs, range0): (_, i32, f64) = arg_utils::parse_3_args(args)?;

    let mut out = String::new();

    let func = "tan(x) / x - 1";
    let poly_i = (1..=num_coeffs).map(|i| i * 2).collect::<Vec<_>>();
    let range = (-range0, range0);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 1, "K", &mut out);

    Ok(out)
}

/// Generates an approximation of `asin(x) - x` in `[-0.00001, 0.50001]` with
/// the odd powers `x^3, x^5, ...`.
///
/// Arguments: `fkind num_coeffs`
pub(super) fn gen_asin_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs) = arg_utils::parse_2_args(args)?;

    let mut out = String::new();

    let func = "asin(x) - x";
    let poly_i = (1..=num_coeffs).map(|i| i * 2 + 1).collect::<Vec<_>>();
    let range = (-0.00001, 0.50001);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 0, "K", &mut out);

    Ok(out)
}

/// Generates an approximation of `atan(x) - x` in `[range_start, range_end]`
/// with the odd powers `x^3, x^5, ...`.
///
/// Arguments: `fkind num_coeffs range_start range_end`
pub(super) fn gen_atan_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs, range_start, range_end) = arg_utils::parse_4_args(args)?;

    let mut out = String::new();

    let func = "atan(x) - x";
    let poly_i = (1..=num_coeffs).map(|i| i * 2 + 1).collect::<Vec<_>>();
    let range = (range_start, range_end);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 0, "K", &mut out);

    Ok(out)
}

/// Generates an approximation of `asinh(x) / x - 1` in
/// `[range_start, range_end]` with the even powers `x^2, x^4, ...`.
///
/// Arguments: `fkind num_coeffs range_start range_end`
pub(super) fn gen_asinh_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs, range_start, range_end) = arg_utils::parse_4_args(args)?;

    let mut out = String::new();

    let func = "asinh(x) / x - 1";
    let poly_i = (1..=num_coeffs).map(|i| i * 2).collect::<Vec<_>>();
    let range = (range_start, range_end);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 1, "K", &mut out);

    Ok(out)
}

/// Generates an approximation of `Γ(x + offset)` in `[range_start, range_end]`
/// with the powers `1, x, ..., x^poly_deg`, minimizing the relative error.
///
/// Arguments: `fkind poly_deg offset range_start range_end`
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

/// Generates an approximation of `ln(Γ(x + offset)) / x` in
/// `[range_start, range_end]` with the powers `1, x, ..., x^(poly_deg - 1)`,
/// minimizing the relative error (i.e., `ln(Γ(x + offset))` is approximated
/// with the powers `x, x^2, ..., x^poly_deg`).
///
/// Arguments: `fkind poly_deg offset range_start range_end`
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

/// Generates an approximation of
/// `Γ(1 / x) / ((1 / x + g - 0.5)^(1 / x - 0.5) * exp(-(1 / x + g - 0.5)))`
/// in `[range_start, range_end]` with the powers `1, x, ..., x^poly_deg`.
///
/// This is a Lanczos-like approximation:
/// `Γ(z) ~= (z + g - 0.5)^(z - 0.5) * exp(-(z + g - 0.5)) * P(1 / z)`.
///
/// Arguments: `fkind poly_deg g range_start range_end`
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
