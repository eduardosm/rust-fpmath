use std::fmt::Write as _;

use crate::RunError;

mod arg_utils;
mod consts;
mod data;
mod julia;
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
        "cbrt::consts" => data::cbrt::gen_consts(&args),
        "cbrt::inv_cbrt_poly" => data::cbrt::gen_inv_cbrt_poly(&args),
        "exp::consts" => data::exp::gen_consts(&args),
        "exp::exp_special_poly" => data::exp::gen_exp_special_poly(&args),
        "exp::exp_m1_special_poly" => data::exp::gen_exp_m1_special_poly(&args),
        "exp2::consts" => data::exp2::gen_consts(&args),
        "exp10::consts" => data::exp10::gen_consts(&args),
        "ln::consts" => data::ln::gen_consts(&args),
        "ln::ln_special_poly" => data::ln::gen_ln_special_poly(&args),
        "ln::ln_special_poly_ex" => data::ln::gen_ln_special_poly_ex(&args),
        "log2::consts" => data::log2::gen_consts(&args),
        "log10::consts" => data::log10::gen_consts(&args),
        "gamma::consts" => data::gamma::gen_consts(&args),
        "gamma::ln_gamma_poly" => data::gamma::gen_ln_gamma_poly(&args),
        "gamma::special_poly" => data::gamma::gen_special_poly(&args),
        "reduce_pi_2::consts" => data::reduce_pi_2::gen_consts(&args),
        "reduce_pi_2_large::frac_2_pi_large" => data::reduce_pi_2_large::gen_frac_2_pi_large(&args),
        "reduce_pi_2_large::frac_pi_2_medium" => {
            data::reduce_pi_2_large::gen_frac_pi_2_medium(&args)
        }
        "reduce_90_deg::consts" => data::reduce_90_deg::gen_consts(&args),
        "reduce_half_mul_pi::consts" => data::reduce_half_mul_pi::gen_consts(&args),
        "sin_cos::consts" => data::sin_cos::gen_consts(&args),
        "sin_cos::sin_poly" => data::sin_cos::gen_sin_poly(&args),
        "sin_cos::sin_poly_ex" => data::sin_cos::gen_sin_poly_ex(&args),
        "sin_cos::cos_poly" => data::sin_cos::gen_cos_poly(&args),
        "tan::tan_poly" => data::tan::gen_tan_poly(&args),
        "rad_to_deg::consts" => data::rad_to_deg::gen_consts(&args),
        "div_pi::consts" => data::div_pi::gen_consts(&args),
        "asin_acos::consts" => data::asin_acos::gen_consts(&args),
        "asin_acos::asin_poly" => data::asin_acos::gen_asin_poly(&args),
        "atan::consts" => data::atan::gen_consts(&args),
        "atan::atan_poly" => data::atan::gen_atan_poly(&args),

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
    NormDoubleF32,
    NormDoubleF64,
    SemiDoubleF32,
    SemiDoubleF64,
}

impl std::str::FromStr for FloatKind {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "f32" => Ok(Self::F32),
            "f64" => Ok(Self::F64),
            "NormDouble<f32>" => Ok(Self::NormDoubleF32),
            "NormDouble<f64>" => Ok(Self::NormDoubleF64),
            "SemiDouble<f32>" => Ok(Self::SemiDoubleF32),
            "SemiDouble<f64>" => Ok(Self::SemiDoubleF64),
            _ => Err("invalid float kind"),
        }
    }
}

impl FloatKind {
    fn rug_aux_prec(self) -> u32 {
        match self {
            Self::F32 => 128,
            Self::F64 => 256,
            Self::NormDoubleF32 => 256,
            Self::NormDoubleF64 => 384,
            Self::SemiDoubleF32 => 192,
            Self::SemiDoubleF64 => 320,
        }
    }

    fn float_prec(self) -> u32 {
        match self {
            Self::F32 => 24,
            Self::F64 => 53,
            Self::NormDoubleF32 => 24 * 2,
            Self::NormDoubleF64 => 53 * 2,
            Self::SemiDoubleF32 => 12 + 24,
            Self::SemiDoubleF64 => 26 + 53,
        }
    }

    fn split_prec(self) -> u32 {
        match self {
            Self::F32 => 12,
            Self::F64 => 26,
            Self::NormDoubleF32 => unimplemented!(),
            Self::NormDoubleF64 => unimplemented!(),
            Self::SemiDoubleF32 => unimplemented!(),
            Self::SemiDoubleF64 => unimplemented!(),
        }
    }

    fn to_norm_double(self) -> Self {
        match self {
            Self::F32 => Self::NormDoubleF32,
            Self::F64 => Self::NormDoubleF64,
            Self::NormDoubleF32 => unimplemented!(),
            Self::NormDoubleF64 => unimplemented!(),
            Self::SemiDoubleF32 => unimplemented!(),
            Self::SemiDoubleF64 => unimplemented!(),
        }
    }

    fn to_semi_double(self) -> Self {
        match self {
            Self::F32 => Self::SemiDoubleF32,
            Self::F64 => Self::SemiDoubleF64,
            Self::NormDoubleF32 => unimplemented!(),
            Self::NormDoubleF64 => unimplemented!(),
            Self::SemiDoubleF32 => unimplemented!(),
            Self::SemiDoubleF64 => unimplemented!(),
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::F32 => "f32",
            Self::F64 => "f64",
            Self::NormDoubleF32 => "NormDouble<f32>",
            Self::NormDoubleF64 => "NormDouble<f64>",
            Self::SemiDoubleF32 => "SemiDouble<f32>",
            Self::SemiDoubleF64 => "SemiDouble<f64>",
        }
    }
}

fn split_hi_lo(mut tmp: rug::Float, hi_prec: u32) -> (rug::Float, rug::Float) {
    let (hi, _) = rug::Float::with_val_round(hi_prec, &tmp, rug::float::Round::Zero);
    tmp -= &hi;
    (hi, tmp)
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
        FloatKind::NormDoubleF32 => {
            let hi = val.to_f32_round(rug::float::Round::Zero);
            let tmp = rug::Float::with_val(val.prec(), val - hi);
            let lo = tmp.to_f32();
            write!(
                out,
                "NormDouble::with_parts(f32::from_bits(0x{:08X}), f32::from_bits(0x{:08X}))",
                hi.to_bits(),
                lo.to_bits(),
            )
            .unwrap();
        }
        FloatKind::NormDoubleF64 => {
            let hi = val.to_f64_round(rug::float::Round::Zero);
            let tmp = rug::Float::with_val(val.prec(), val - hi);
            let lo = tmp.to_f64();
            write!(
                out,
                "NormDouble::with_parts(f64::from_bits(0x{:016X}), f64::from_bits(0x{:016X}))",
                hi.to_bits(),
                lo.to_bits(),
            )
            .unwrap();
        }
        FloatKind::SemiDoubleF32 => {
            let (hi, lo) = split_hi_lo(val.clone(), 12);
            let hi = hi.to_f32_round(rug::float::Round::Zero);
            let lo = lo.to_f32();
            write!(
                out,
                "SemiDouble::with_parts(f32::from_bits(0x{:08X}), f32::from_bits(0x{:08X}))",
                hi.to_bits(),
                lo.to_bits(),
            )
            .unwrap();
        }
        FloatKind::SemiDoubleF64 => {
            let (hi, lo) = split_hi_lo(val.clone(), 26);
            let hi = hi.to_f64_round(rug::float::Round::Zero);
            let lo = lo.to_f64();
            write!(
                out,
                "SemiDouble::with_parts(f64::from_bits(0x{:016X}), f64::from_bits(0x{:016X}))",
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
        FloatKind::NormDoubleF32
        | FloatKind::NormDoubleF64
        | FloatKind::SemiDoubleF32
        | FloatKind::SemiDoubleF64 => {
            let prec = fkind.float_prec();
            let dec_prec = ((prec as f64) * std::f64::consts::LOG10_2) as usize;
            write!(out, "{val:.dec_prec$e}").unwrap();
        }
    }
}
