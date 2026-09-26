use std::fmt::Write as _;

use crate::RunError;

mod approx;
mod arg_utils;
mod consts;
mod julia;
mod reduce_pi_2;
mod sollya;

pub(crate) fn generate(param: &str) -> Result<String, RunError> {
    let mut args = param.split_ascii_whitespace();
    let cmd = args.next().ok_or_else(|| {
        eprintln!("Empty generate parameter");
        RunError
    })?;
    let args = args.collect::<Vec<_>>();

    let r = match cmd {
        "consts" => consts::gen_consts(&args),

        "reduce_pi_2::medium_consts" => reduce_pi_2::gen_medium_consts(&args),
        "reduce_pi_2::frac_2_pi_large" => reduce_pi_2::gen_frac_2_pi_large(&args),
        "reduce_pi_2::frac_pi_2_medium" => reduce_pi_2::gen_frac_pi_2_medium(&args),

        "inv_cbrt_poly" => approx::gen_inv_cbrt_poly(&args),
        "exp_m1_poly" => approx::gen_exp_m1_poly(&args),
        "ln_1p_poly" => approx::gen_ln_1p_poly(&args),
        "ln_table" => approx::gen_ln_table(&args),
        "ln_lo_scale_table" => approx::gen_ln_lo_scale_table(&args),
        "sin_poly" => approx::gen_sin_poly(&args),
        "cos_poly" => approx::gen_cos_poly(&args),
        "tan_poly" => approx::gen_tan_poly(&args),
        "asin_poly" => approx::gen_asin_poly(&args),
        "atan_poly" => approx::gen_atan_poly(&args),
        "asinh_poly" => approx::gen_asinh_poly(&args),
        "gamma_poly" => approx::gen_gamma_poly(&args),
        "ln_gamma_poly" => approx::gen_ln_gamma_poly(&args),
        "gamma_lanczos_poly" => approx::gen_gamma_lanczos_poly(&args),

        _ => {
            eprintln!("Invalid generate parameter: {cmd:?}");
            return Err(RunError);
        }
    };

    r.map_err(|e| {
        eprintln!("Error generating {cmd:?}: {e}");
        RunError
    })
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum FloatKind {
    F32,
    F64,
    F64x2,
}

impl std::str::FromStr for FloatKind {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "f32" => Ok(Self::F32),
            "f64" => Ok(Self::F64),
            "F64x2" => Ok(Self::F64x2),
            _ => Err("invalid float kind"),
        }
    }
}

impl FloatKind {
    fn rug_aux_prec(self) -> u32 {
        match self {
            Self::F32 => 128,
            Self::F64 => 256,
            Self::F64x2 => 384,
        }
    }

    fn float_prec(self) -> u32 {
        match self {
            Self::F32 => 24,
            Self::F64 => 53,
            Self::F64x2 => 53 * 2,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::F32 => "f32",
            Self::F64 => "f64",
            Self::F64x2 => "F64x2",
        }
    }
}

fn render_const(fkind: FloatKind, name: &str, val: rug::Float, out: &mut String) {
    out.push_str("const ");
    out.push_str(name);
    out.push_str(": ");
    out.push_str(fkind.name());
    out.push_str(" = ");
    render_const_value(fkind, &val, out);
    out.push_str("; // ");
    render_const_dec_value(fkind, &val, out);
    out.push('\n');
}

fn render_const_value(fkind: FloatKind, val: &rug::Float, out: &mut String) {
    match fkind {
        FloatKind::F32 => {
            let val = val.to_f32();
            write!(out, "f32::from_bits(0x{:08X})", val.to_bits()).unwrap();
        }
        FloatKind::F64 => {
            let val = val.to_f64();
            write!(out, "f64::from_bits(0x{:016X})", val.to_bits()).unwrap();
        }
        FloatKind::F64x2 => {
            let hi = val.to_f64_round(rug::float::Round::Nearest);
            let tmp = rug::Float::with_val(val.prec(), val - hi);
            let lo = tmp.to_f64();
            write!(
                out,
                "F64x2::from_bits(0x{:016X}, 0x{:016X})",
                hi.to_bits(),
                lo.to_bits(),
            )
            .unwrap();
        }
    }
}

fn render_const_dec_value(fkind: FloatKind, val: &rug::Float, out: &mut String) {
    match fkind {
        FloatKind::F32 => {
            let val = val.to_f32();
            write!(out, "{val:e}").unwrap();
        }
        FloatKind::F64 => {
            let val = val.to_f64();
            write!(out, "{val:e}").unwrap();
        }
        FloatKind::F64x2 => {
            let prec = fkind.float_prec();
            let dec_prec = ((prec as f64) * std::f64::consts::LOG10_2) as usize;
            write!(out, "{val:.dec_prec$e}").unwrap();
        }
    }
}
