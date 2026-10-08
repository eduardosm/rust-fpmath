use std::fmt::Write as _;

use super::{
    CoeffKinds, FloatKind, arg_utils, julia, render_const, render_const_dec_value,
    render_const_value, sollya,
};

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
    let (args, num_fixed) = split_num_fixed(args)?;
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

/// Splits the optional last argument `num_fixed` (number of fixed leading
/// coefficients of a series, at least 1, and 1 when not present) from the
/// four arguments `fkind num_coeffs range_start range_end`.
fn split_num_fixed<'a, 'b>(args: &'a [&'b str]) -> Result<(&'a [&'b str], i32), String> {
    match args {
        [args @ .., num_fixed] if args.len() == 4 => {
            let num_fixed: i32 = num_fixed
                .parse()
                .map_err(|e| format!("failed to parse fifth argument {num_fixed:?}: {e}"))?;
            if num_fixed < 1 {
                return Err(format!("invalid number of fixed terms: {num_fixed}"));
            }
            Ok((args, num_fixed))
        }
        [_, _, _, _] => Ok((args, 1)),
        _ => Err(format!("expected 4 or 5 arguments, found {}", args.len())),
    }
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
/// When `num_fixed` (default 1) is greater than 1, the coefficients of
/// `x^1, ..., x^(num_fixed - 1)` are fixed to those of the Taylor series
/// (`-1/2, 1/3, ...`), so it approximates
/// `ln(1 + x) / x - 1 + x / 2 - ... + (-1)^num_fixed * x^(num_fixed - 1) / num_fixed`
/// with the powers `x^num_fixed, x^(num_fixed + 1), ...`.
///
/// Arguments: `fkind num_coeffs range_start range_end [num_fixed]`
pub(super) fn gen_ln_1p_poly(args: &[&str]) -> Result<String, String> {
    let (args, num_fixed) = split_num_fixed(args)?;
    let (fkind, num_coeffs, range_start, range_end): (_, i32, _, _) =
        arg_utils::parse_4_args(args)?;

    let mut out = String::new();

    let mut func = String::from("log1p(x) / x - 1");
    for i in 2..=num_fixed {
        let sign = if i % 2 == 0 { '+' } else { '-' };
        write!(func, " {sign} x^{} / {i}", i - 1).unwrap();
    }
    let poly_i = (num_fixed..(num_fixed + num_coeffs)).collect::<Vec<_>>();
    let range = (range_start, range_end);

    sollya::run_and_render_remez(fkind, &func, range, &poly_i, 1, "K", &mut out);

    Ok(out)
}

/// Generates a table of triples `(s, hi, lo)` of raw `f64` bits, for `i` in
/// `0..=2^bits`, where `s` is `1 / (1 + i / 2^bits)` rounded to
/// `scale_fkind`, and `hi + lo = -ln(s)` (approximately `ln(1 + i / 2^bits)`),
/// with `hi` rounded to a multiple of `2^-hi_exp` and `lo` being the
/// remainder rounded to `f64`.
///
/// Arguments: `scale_fkind bits hi_exp`
pub(super) fn gen_ln_table(args: &[&str]) -> Result<String, String> {
    let (scale_fkind, bits, hi_exp): (FloatKind, u32, u32) = arg_utils::parse_3_args(args)?;

    let mut out = String::new();

    let prec = FloatKind::F64x2.rug_aux_prec();
    let scale_prec = scale_fkind.float_prec();
    let scale_ftype = scale_fkind.name();
    let num = 1 << bits;

    writeln!(
        out,
        "// LN_TBL[i] = (bits(s), bits(hi), bits(lo)), where s = 1 / (1 + i / {num}) rounded to {scale_ftype}",
    )
    .unwrap();
    writeln!(
        out,
        "// and hi + lo = -ln(s), with hi being a multiple of 2^-{hi_exp}",
    )
    .unwrap();
    writeln!(out, "static LN_TBL: [(u64, u64, u64); {}] = [", num + 1).unwrap();
    for x in 0..=num {
        let v = ((rug::Float::with_val(prec, x) >> bits) + 1u8).recip();
        // The logarithm is calculated from the rounded reciprocal, which is
        // the one used in the reduction.
        let (v, _) = rug::Float::with_val_round(scale_prec, v, rug::float::Round::Nearest);
        let v = rug::Float::with_val(prec, v);
        let ln_v = -v.clone().ln();
        let hi = rug::Float::with_val(prec, &ln_v << hi_exp).round() >> hi_exp;
        let hi_f64 = hi.to_f64();
        assert!(hi == hi_f64, "{hi} is not exactly representable as f64");
        let lo_f64 = rug::Float::with_val(prec, &ln_v - &hi).to_f64();
        write!(
            out,
            "    (0x{:016X}, 0x{:016X}, 0x{:016X}), // ",
            v.to_f64().to_bits(),
            hi_f64.to_bits(),
            lo_f64.to_bits(),
        )
        .unwrap();
        render_const_dec_value(FloatKind::F64x2, &ln_v, &mut out);
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

/// Generates an approximation of `atanh(x) / x - 1` in
/// `[range_start, range_end]` with the even powers `x^2, x^4, ...`.
///
/// Arguments: `fkind num_coeffs range_start range_end`
pub(super) fn gen_atanh_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, num_coeffs, range_start, range_end) = arg_utils::parse_4_args(args)?;

    let mut out = String::new();

    let func = "atanh(x) / x - 1";
    let poly_i = (1..=num_coeffs).map(|i| i * 2).collect::<Vec<_>>();
    let range = (range_start, range_end);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 1, "K", &mut out);

    Ok(out)
}

/// Generates an approximation of `(g(x) - s * x / 4) / x^2` in
/// `[range_start, range_end]` with the powers `1, x, x^2, ...`, where
/// `g(x) = ln((1 + sqrt(1 + s * x)) / 2)`, with `s = 1` for `asinh` and
/// `s = -1` for `acosh`. The coefficients are named after the powers of
/// `g(x) = s * x / 4 + K2 * x^2 + K3 * x^3 + ...`.
///
/// For large `y`, `asinh(y) = ln(2 * y) + g(1 / y^2)` (with `s = 1`) and
/// `acosh(y) = ln(2 * y) + g(1 / y^2)` (with `s = -1`).
///
/// The approximated expression is not defined at 0, so `range_start` must be
/// positive (a tiny value, such as `1e-30`, is equivalent to 0).
///
/// Arguments: `fkind (asinh|acosh) num_coeffs range_start range_end`
pub(super) fn gen_asinh_acosh_large_poly(args: &[&str]) -> Result<String, String> {
    let (fkind, func, num_coeffs, range_start, range_end): (_, String, i32, _, _) =
        arg_utils::parse_5_args(args)?;

    let func = match func.as_str() {
        "asinh" => "(log((1 + sqrt(1 + x)) / 2) - x / 4) / x^2",
        "acosh" => "(log((1 + sqrt(1 - x)) / 2) + x / 4) / x^2",
        _ => return Err(format!("invalid function: {func:?}")),
    };

    let mut out = String::new();

    let poly_i = (0..num_coeffs).collect::<Vec<_>>();
    let range = (range_start, range_end);

    sollya::run_and_render_remez(fkind, func, range, &poly_i, 2, "K", &mut out);

    Ok(out)
}

/// Generates an approximation of `Γ(x + offset)` in `[range_start, range_end]`
/// with the powers `1, x, ..., x^poly_deg`, minimizing the relative error.
///
/// Arguments: `fkind poly_deg offset range_start range_end`
pub(super) fn gen_gamma_poly(args: &[&str]) -> Result<String, String> {
    let (fkinds, poly_deg, offset, range_start, range_end): (CoeffKinds, i32, f64, f64, f64) =
        arg_utils::parse_5_args(args)?;

    let mut out = String::new();

    let func = format!("SpecialFunctions.gamma(x + {offset})");
    let wfunc = "1 / fx";
    let range = (range_start, range_end);

    julia::run_and_render_remez(&fkinds, &func, wfunc, range, poly_deg, 0, "K", &mut out);

    Ok(out)
}

/// Generates an approximation of `ln(Γ(x + offset)) / x` in
/// `[range_start, range_end]` with the powers `1, x, ..., x^(poly_deg - 1)`,
/// minimizing the relative error (i.e., `ln(Γ(x + offset))` is approximated
/// with the powers `x, x^2, ..., x^poly_deg`).
///
/// Arguments: `fkind poly_deg offset range_start range_end`
pub(super) fn gen_ln_gamma_poly(args: &[&str]) -> Result<String, String> {
    let (fkinds, poly_deg, offset, range_start, range_end): (CoeffKinds, i32, f64, f64, f64) =
        arg_utils::parse_5_args(args)?;

    let mut out = String::new();

    let func = format!("SpecialFunctions.lgamma(x + {offset}) / x");
    let wfunc = "1 / fx";
    let range = (range_start, range_end);

    julia::run_and_render_remez(&fkinds, &func, wfunc, range, poly_deg - 1, 1, "K", &mut out);

    Ok(out)
}

/// Generates an expansion of `ln(|Γ(x)|)` around one of its negative roots
/// `x0`, as a constant of type `LnGammaRoot`:
///
/// `ln(|Γ(x0 + h)|) ~= h * (k0 + h * (k[0] + h * (k[1] + ...)))`
///
/// fitted in `[-2^radius_exp, 2^radius_exp]` minimizing the relative error,
/// with `x0` split in three `f64` (`x0[0] + x0[1] + x0[2]`), `k0` as `F64x2`
/// and the other coefficients as an array of `num_coeffs - 1` `f64`.
///
/// Arguments: `name x0_approx radius_exp num_coeffs`, where `x0_approx` is an
/// approximation of the root.
pub(super) fn gen_ln_gamma_root(args: &[&str]) -> Result<String, String> {
    let (name, x0_approx, radius_exp, num_coeffs): (String, f64, i32, i32) =
        arg_utils::parse_4_args(args)?;
    if num_coeffs < 2 {
        return Err("expected at least 2 coefficients".into());
    }

    // Refine the root with Newton's method
    let root_prec = 2048;
    let mut x0 = rug::Float::with_val(root_prec, x0_approx);
    for _ in 0..20 {
        let step = x0.clone().ln_abs_gamma().0 / x0.clone().digamma();
        x0 -= step;
    }
    let ln_gamma_x0 = x0.clone().ln_abs_gamma().0;
    if !ln_gamma_x0.is_zero() && ln_gamma_x0.get_exp().unwrap() > -(root_prec as i32 - 64) {
        return Err(format!("root near {x0_approx} did not converge"));
    }
    if (x0.clone() - x0_approx).abs() > 1e-12 {
        return Err(format!("root near {x0_approx} converged to {x0:e}"));
    }

    // The expansion is used in [x0 - r, x0 + r], which must not contain a pole
    // and must be within a single binade, so subtracting `x0[0]` is exact.
    let radius = rug::Float::with_val(64, rug::Float::i_exp(1, radius_exp));
    let range_start = x0.clone() - &radius;
    let range_end = x0.clone() + &radius;
    if range_end >= range_start.clone().ceil() {
        return Err(format!("range around {x0_approx} contains a pole"));
    }
    if range_start.get_exp() != range_end.get_exp() {
        return Err(format!("range around {x0_approx} spans multiple binades"));
    }

    // x0 = x0_0 + x0_1 + x0_2
    let mut x0_parts = Vec::new();
    let mut rem = x0.clone();
    for _ in 0..3 {
        let part = rem.to_f64();
        x0_parts.push(part);
        rem -= part;
    }

    // Fit ln(|Γ(x0 + h)|) / h, minimizing the relative error
    let x0_digits = 1024 * 3 / 10 + 10;
    let x0_str = x0.to_string_radix(10, Some(x0_digits));
    let x0_big = format!("BigFloat(\"{x0_str}\")");
    let func = format!(
        "iszero(x) ? SpecialFunctions.digamma({x0_big}) : \
        SpecialFunctions.logabsgamma({x0_big} + x)[1] / x"
    );
    let wfunc = "1 / abs(fx)";
    let r = radius.to_f64();
    let coeffs = julia::run_remez(&func, wfunc, (-r, r), num_coeffs - 1);

    let mut out = String::new();
    let render_value = |fkind: FloatKind, value: &rug::Float, out: &mut String| {
        render_const_value(fkind, value, out);
        out.push_str(", // ");
        render_const_dec_value(fkind, value, out);
        out.push('\n');
    };

    writeln!(out, "const {name}: LnGammaRoot = LnGammaRoot {{").unwrap();
    writeln!(out, "x0: [").unwrap();
    for part in x0_parts {
        render_value(FloatKind::F64, &rug::Float::with_val(53, part), &mut out);
    }
    writeln!(out, "],").unwrap();
    out.push_str("radius: ");
    render_value(FloatKind::F64, &radius, &mut out);
    out.push_str("k0: ");
    render_value(FloatKind::F64x2, &coeffs[0], &mut out);
    writeln!(out, "k: [").unwrap();
    for coeff in coeffs[1..].iter() {
        render_value(FloatKind::F64, coeff, &mut out);
    }
    writeln!(out, "],").unwrap();
    writeln!(out, "}};").unwrap();

    Ok(out)
}

/// Generates the coefficients of the Stirling series of `ln(Γ(x))`, for
/// `x >= x0`:
///
/// `ln(Γ(x)) ~= (x - 0.5) * ln(x) - x + K0 + K1 / x + K3 / x^3 + K5 / x^5 + ...`
///
/// where `K0 = ln(2 * π) / 2` and `K1 = 1/12` (as in the series), and
/// `K3, K5, ...` are fitted to minimize the absolute error.
///
/// Arguments: `fkinds num_coeffs x0`, where `num_coeffs` includes `K0` and
/// `K1`.
pub(super) fn gen_gamma_stirling_poly(args: &[&str]) -> Result<String, String> {
    let (fkinds, num_coeffs, x0): (CoeffKinds, i32, f64) = arg_utils::parse_3_args(args)?;
    if num_coeffs < 3 {
        return Err("expected at least 3 coefficients".into());
    }

    let mut out = String::new();

    // With t = 1 / x and s = t^2, fit
    // R(s) = (ln(Γ(1 / t)) - (1 / t - 0.5) * ln(1 / t) + 1 / t - K0 - K1 * t) / t^3
    // with a polynomial in `s`, minimizing the absolute error of `t^3 * R(s)`.
    let func = "begin \
        t = sqrt(x); \
        (SpecialFunctions.loggamma(1 / t) - (1 / t - 0.5) * log(1 / t) + 1 / t \
        - log(2 * BigFloat(pi)) / 2 - t / 12) / (t * t * t) \
        end";
    let wfunc = "x * sqrt(x)";
    let range = (1e-40, 1.0 / (x0 * x0));

    let coeffs = julia::run_remez(func, wfunc, range, num_coeffs - 3);

    let prec = 1024;
    let pi = rug::Float::with_val(prec, rug::float::Constant::Pi);
    let k0 = (pi * 2u8).ln() / 2u8;
    render_const(fkinds.get(0), "K0", k0, &mut out);
    let k1 = rug::Float::with_val(prec, 1) / 12u8;
    render_const(fkinds.get(1), "K1", k1, &mut out);
    for (j, coeff) in coeffs.into_iter().enumerate() {
        render_const(
            fkinds.get(j + 2),
            &format!("K{}", j * 2 + 3),
            coeff,
            &mut out,
        );
    }

    Ok(out)
}
