use crate::RunError;

mod common;
mod f32;
mod f64;
mod julia;
mod sollya;

pub(crate) fn generate(param: &str) -> Result<String, RunError> {
    let out = match param {
        // common
        "common::reduce_pi_2_large::FRAC_2_PI_LARGE" => {
            common::reduce_pi_2_large::gen_frac_2_pi_large()
        }
        "common::reduce_pi_2_large::FRAC_PI_2_MEDIUM" => {
            common::reduce_pi_2_large::gen_frac_pi_2_medium()
        }

        // f32
        "f32::consts" => f32::gen_consts(),
        "f32::cbrt::consts" => f32::cbrt::gen_consts(),
        "f32::cbrt::inv_cbrt_poly" => f32::cbrt::gen_inv_cbrt_poly(),
        "f32::div_pi::consts" => f32::div_pi::gen_consts(),
        "f32::exp::consts" => f32::exp::gen_consts(),
        "f32::exp::exp_special_poly" => f32::exp::gen_exp_special_poly(),
        "f32::exp::exp_m1_special_poly" => f32::exp::gen_exp_m1_special_poly(),
        "f32::exp2::consts" => f32::exp2::gen_consts(),
        "f32::exp10::consts" => f32::exp10::gen_consts(),
        "f32::gamma::consts" => f32::gamma::gen_consts(),
        "f32::gamma::lgamma_poly_1" => f32::gamma::gen_lgamma_poly_1(),
        "f32::gamma::lgamma_poly_2" => f32::gamma::gen_lgamma_poly_2(),
        "f32::gamma::special_poly" => f32::gamma::gen_special_poly(),
        "f32::log::consts" => f32::log::gen_consts(),
        "f32::log::log_special_poly" => f32::log::gen_log_special_poly(),
        "f32::log::log_special_poly_ex" => f32::log::gen_log_special_poly_ex(),
        "f32::log2::consts" => f32::log2::gen_consts(),
        "f32::log10::consts" => f32::log10::gen_consts(),
        "f32::rad_to_deg::consts" => f32::rad_to_deg::gen_consts(),
        "f32::reduce_90_deg::consts" => f32::reduce_90_deg::gen_consts(),
        "f32::reduce_half_mul_pi::consts" => f32::reduce_half_mul_pi::gen_consts(),
        "f32::reduce_pi_2::consts" => f32::reduce_pi_2::gen_consts(),
        "f32::sin_cos::consts" => f32::sin_cos::gen_consts(),
        "f32::sin_cos::sin_poly" => f32::sin_cos::gen_sin_poly(),
        "f32::sin_cos::sin_poly_ex" => f32::sin_cos::gen_sin_poly_ex(),
        "f32::sin_cos::cos_poly" => f32::sin_cos::gen_cos_poly(),
        "f32::tan::tan_poly" => f32::tan::gen_tan_poly(),
        "f32::asin_acos::consts" => f32::asin_acos::gen_consts(),
        "f32::asin_acos::asin_poly" => f32::asin_acos::gen_asin_poly(),
        "f32::atan::consts" => f32::atan::gen_consts(),
        "f32::atan::atan_poly" => f32::atan::gen_atan_poly(),

        // f64
        "f64::consts" => f64::gen_consts(),
        "f64::cbrt::consts" => f64::cbrt::gen_consts(),
        "f64::cbrt::inv_cbrt_poly" => f64::cbrt::gen_inv_cbrt_poly(),
        "f64::div_pi::consts" => f64::div_pi::gen_consts(),
        "f64::exp::consts" => f64::exp::gen_consts(),
        "f64::exp::exp_special_poly" => f64::exp::gen_exp_special_poly(),
        "f64::exp::exp_m1_special_poly" => f64::exp::gen_exp_m1_special_poly(),
        "f64::exp2::consts" => f64::exp2::gen_consts(),
        "f64::exp10::consts" => f64::exp10::gen_consts(),
        "f64::gamma::consts" => f64::gamma::gen_consts(),
        "f64::gamma::lgamma_poly_1" => f64::gamma::gen_lgamma_poly_1(),
        "f64::gamma::lgamma_poly_2" => f64::gamma::gen_lgamma_poly_2(),
        "f64::gamma::special_poly" => f64::gamma::gen_special_poly(),
        "f64::log::consts" => f64::log::gen_consts(),
        "f64::log::log_special_poly" => f64::log::gen_log_special_poly(),
        "f64::log::log_special_poly_ex" => f64::log::gen_log_special_poly_ex(),
        "f64::log2::consts" => f64::log2::gen_consts(),
        "f64::log10::consts" => f64::log10::gen_consts(),
        "f64::rad_to_deg::consts" => f64::rad_to_deg::gen_consts(),
        "f64::reduce_90_deg::consts" => f64::reduce_90_deg::gen_consts(),
        "f64::reduce_half_mul_pi::consts" => f64::reduce_half_mul_pi::gen_consts(),
        "f64::reduce_pi_2::consts" => f64::reduce_pi_2::gen_consts(),
        "f64::sin_cos::consts" => f64::sin_cos::gen_consts(),
        "f64::sin_cos::sin_poly" => f64::sin_cos::gen_sin_poly(),
        "f64::sin_cos::sin_poly_ex" => f64::sin_cos::gen_sin_poly_ex(),
        "f64::sin_cos::cos_poly" => f64::sin_cos::gen_cos_poly(),
        "f64::tan::tan_poly" => f64::tan::gen_tan_poly(),
        "f64::asin_acos::consts" => f64::asin_acos::gen_consts(),
        "f64::asin_acos::asin_poly" => f64::asin_acos::gen_asin_poly(),
        "f64::atan::consts" => f64::atan::gen_consts(),
        "f64::atan::atan_poly" => f64::atan::gen_atan_poly(),

        _ => {
            eprintln!("Invalid generate parameter: {param:?}");
            return Err(RunError);
        }
    };
    Ok(out)
}
