use std::fmt::Write as _;

use crate::RunError;

mod arg_utils;
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
        "consts" => data::gen_consts(&args),
        "cbrt::consts" => data::cbrt::gen_consts(&args),
        "cbrt::inv_cbrt_poly" => data::cbrt::gen_inv_cbrt_poly(&args),
        "exp::consts" => data::exp::gen_consts(&args),
        "exp::exp_special_poly" => data::exp::gen_exp_special_poly(&args),
        "exp::exp_m1_special_poly" => data::exp::gen_exp_m1_special_poly(&args),
        "exp2::consts" => data::exp2::gen_consts(&args),
        "exp10::consts" => data::exp10::gen_consts(&args),
        "log::consts" => data::log::gen_consts(&args),
        "log::log_special_poly" => data::log::gen_log_special_poly(&args),
        "log::log_special_poly_ex" => data::log::gen_log_special_poly_ex(&args),
        "log2::consts" => data::log2::gen_consts(&args),
        "log10::consts" => data::log10::gen_consts(&args),
        "gamma::consts" => data::gamma::gen_consts(&args),
        "gamma::lgamma_poly" => data::gamma::gen_lgamma_poly(&args),
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
}

impl std::str::FromStr for FloatKind {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "f32" => Ok(Self::F32),
            "f64" => Ok(Self::F64),
            _ => Err("invalid float kind"),
        }
    }
}

impl FloatKind {
    fn rug_aux_prec(self) -> u32 {
        match self {
            Self::F32 => 128,
            Self::F64 => 256,
        }
    }

    fn float_prec(self) -> u32 {
        match self {
            Self::F32 => 24,
            Self::F64 => 53,
        }
    }

    fn split_prec(self) -> u32 {
        match self {
            Self::F32 => 12,
            Self::F64 => 26,
        }
    }
}

fn split_hi_lo(mut tmp: rug::Float, hi_prec: u32) -> (rug::Float, rug::Float) {
    let (hi, _) = rug::Float::with_val_round(hi_prec, &tmp, rug::float::Round::Zero);
    tmp -= &hi;
    (hi, tmp)
}

fn render_const(fkind: FloatKind, name: &str, val: rug::Float, out: &mut String) {
    match fkind {
        FloatKind::F32 => {
            let val = val.to_f32();
            writeln!(
                out,
                "const {name}: u32 = 0x{:08X}; // {val:e}",
                val.to_bits(),
            )
            .unwrap();
        }
        FloatKind::F64 => {
            let val = val.to_f64();
            writeln!(
                out,
                "const {name}: u64 = 0x{:016X}; // {val:e}",
                val.to_bits(),
            )
            .unwrap();
        }
    }
}
