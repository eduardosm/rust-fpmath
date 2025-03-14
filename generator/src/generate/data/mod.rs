use super::{arg_utils, render_const, FloatKind};

pub(super) mod asin_acos;
pub(super) mod atan;
pub(super) mod cbrt;
pub(super) mod div_pi;
pub(super) mod exp;
pub(super) mod exp10;
pub(super) mod exp2;
pub(super) mod gamma;
pub(super) mod log;
pub(super) mod log10;
pub(super) mod log2;
pub(super) mod rad_to_deg;
pub(super) mod reduce_90_deg;
pub(super) mod reduce_half_mul_pi;
pub(super) mod reduce_pi_2;
pub(super) mod reduce_pi_2_large;
pub(super) mod sin_cos;
pub(super) mod tan;

pub(super) fn gen_consts(args: &[&str]) -> Result<String, String> {
    let fkind: FloatKind = arg_utils::parse_1_arg(args)?;
    let aux_prec = fkind.rug_aux_prec();

    let mut out = String::new();

    // π
    let tmp = rug::Float::with_val(aux_prec, rug::float::Constant::Pi);
    render_const(fkind, "PI", tmp, &mut out);

    // π/2
    let tmp = rug::Float::with_val(aux_prec, rug::float::Constant::Pi) / 2;
    render_const(fkind, "FRAC_PI_2", tmp, &mut out);

    // π/4
    let tmp = rug::Float::with_val(aux_prec, rug::float::Constant::Pi) / 4;
    render_const(fkind, "FRAC_PI_4", tmp, &mut out);

    // 2/π
    let tmp = 2u8 / rug::Float::with_val(aux_prec, rug::float::Constant::Pi);
    render_const(fkind, "FRAC_2_PI", tmp, &mut out);

    Ok(out)
}
