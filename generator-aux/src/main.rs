#![warn(
    rust_2018_idioms,
    trivial_casts,
    trivial_numeric_casts,
    unreachable_pub,
    unused_qualifications
)]
#![forbid(unsafe_code)]

use std::process::ExitCode;

mod generate;
mod julia;
mod sollya;

fn main() -> ExitCode {
    let args: Vec<_> = std::env::args_os().collect();
    if args.len() != 2 {
        eprintln!("Usage: {} <param>", args[0].to_string_lossy());
        return ExitCode::FAILURE;
    }

    let param = args[1].as_os_str();
    match param.to_str() {
        // common
        Some("common::reduce_pi_2_large::FRAC_2_PI_LARGE") => {
            generate::common::reduce_pi_2_large::gen_frac_2_pi_large();
        }
        Some("common::reduce_pi_2_large::FRAC_PI_2_MEDIUM") => {
            generate::common::reduce_pi_2_large::gen_frac_pi_2_medium();
        }

        // f32
        Some("f32::consts") => generate::f32::gen_consts(),
        Some("f32::cbrt::consts") => generate::f32::cbrt::gen_consts(),
        Some("f32::cbrt::inv_cbrt_poly") => generate::f32::cbrt::gen_inv_cbrt_poly(),
        Some("f32::div_pi::consts") => generate::f32::div_pi::gen_consts(),
        Some("f32::exp::consts") => generate::f32::exp::gen_consts(),
        Some("f32::exp::exp_special_poly") => generate::f32::exp::gen_exp_special_poly(),
        Some("f32::exp::exp_m1_special_poly") => generate::f32::exp::gen_exp_m1_special_poly(),
        Some("f32::exp2::consts") => generate::f32::exp2::gen_consts(),
        Some("f32::exp10::consts") => generate::f32::exp10::gen_consts(),
        Some("f32::gamma::consts") => generate::f32::gamma::gen_consts(),
        Some("f32::gamma::lgamma_poly_1") => generate::f32::gamma::gen_lgamma_poly_1(),
        Some("f32::gamma::lgamma_poly_2") => generate::f32::gamma::gen_lgamma_poly_2(),
        Some("f32::gamma::special_poly") => generate::f32::gamma::gen_special_poly(),
        Some("f32::log::consts") => generate::f32::log::gen_consts(),
        Some("f32::log::log_special_poly") => generate::f32::log::gen_log_special_poly(),
        Some("f32::log::log_special_poly_ex") => generate::f32::log::gen_log_special_poly_ex(),
        Some("f32::log2::consts") => generate::f32::log2::gen_consts(),
        Some("f32::log10::consts") => generate::f32::log10::gen_consts(),
        Some("f32::rad_to_deg::consts") => generate::f32::rad_to_deg::gen_consts(),
        Some("f32::reduce_90_deg::consts") => generate::f32::reduce_90_deg::gen_consts(),
        Some("f32::reduce_half_mul_pi::consts") => generate::f32::reduce_half_mul_pi::gen_consts(),
        Some("f32::reduce_pi_2::consts") => generate::f32::reduce_pi_2::gen_consts(),
        Some("f32::sin_cos::consts") => generate::f32::sin_cos::gen_consts(),
        Some("f32::sin_cos::sin_poly") => generate::f32::sin_cos::gen_sin_poly(),
        Some("f32::sin_cos::sin_poly_ex") => generate::f32::sin_cos::gen_sin_poly_ex(),
        Some("f32::sin_cos::cos_poly") => generate::f32::sin_cos::gen_cos_poly(),
        Some("f32::tan::tan_poly") => generate::f32::tan::gen_tan_poly(),
        Some("f32::asin_acos::consts") => generate::f32::asin_acos::gen_consts(),
        Some("f32::asin_acos::asin_poly") => generate::f32::asin_acos::gen_asin_poly(),
        Some("f32::atan::consts") => generate::f32::atan::gen_consts(),
        Some("f32::atan::atan_poly") => generate::f32::atan::gen_atan_poly(),

        // f64
        Some("f64::consts") => generate::f64::gen_consts(),
        Some("f64::cbrt::consts") => generate::f64::cbrt::gen_consts(),
        Some("f64::cbrt::inv_cbrt_poly") => generate::f64::cbrt::gen_inv_cbrt_poly(),
        Some("f64::div_pi::consts") => generate::f64::div_pi::gen_consts(),
        Some("f64::exp::consts") => generate::f64::exp::gen_consts(),
        Some("f64::exp::exp_special_poly") => generate::f64::exp::gen_exp_special_poly(),
        Some("f64::exp::exp_m1_special_poly") => generate::f64::exp::gen_exp_m1_special_poly(),
        Some("f64::exp2::consts") => generate::f64::exp2::gen_consts(),
        Some("f64::exp10::consts") => generate::f64::exp10::gen_consts(),
        Some("f64::gamma::consts") => generate::f64::gamma::gen_consts(),
        Some("f64::gamma::lgamma_poly_1") => generate::f64::gamma::gen_lgamma_poly_1(),
        Some("f64::gamma::lgamma_poly_2") => generate::f64::gamma::gen_lgamma_poly_2(),
        Some("f64::gamma::special_poly") => generate::f64::gamma::gen_special_poly(),
        Some("f64::log::consts") => generate::f64::log::gen_consts(),
        Some("f64::log::log_special_poly") => generate::f64::log::gen_log_special_poly(),
        Some("f64::log::log_special_poly_ex") => generate::f64::log::gen_log_special_poly_ex(),
        Some("f64::log2::consts") => generate::f64::log2::gen_consts(),
        Some("f64::log10::consts") => generate::f64::log10::gen_consts(),
        Some("f64::rad_to_deg::consts") => generate::f64::rad_to_deg::gen_consts(),
        Some("f64::reduce_90_deg::consts") => generate::f64::reduce_90_deg::gen_consts(),
        Some("f64::reduce_half_mul_pi::consts") => generate::f64::reduce_half_mul_pi::gen_consts(),
        Some("f64::reduce_pi_2::consts") => generate::f64::reduce_pi_2::gen_consts(),
        Some("f64::sin_cos::consts") => generate::f64::sin_cos::gen_consts(),
        Some("f64::sin_cos::sin_poly") => generate::f64::sin_cos::gen_sin_poly(),
        Some("f64::sin_cos::sin_poly_ex") => generate::f64::sin_cos::gen_sin_poly_ex(),
        Some("f64::sin_cos::cos_poly") => generate::f64::sin_cos::gen_cos_poly(),
        Some("f64::tan::tan_poly") => generate::f64::tan::gen_tan_poly(),
        Some("f64::asin_acos::consts") => generate::f64::asin_acos::gen_consts(),
        Some("f64::asin_acos::asin_poly") => generate::f64::asin_acos::gen_asin_poly(),
        Some("f64::atan::consts") => generate::f64::atan::gen_consts(),
        Some("f64::atan::atan_poly") => generate::f64::atan::gen_atan_poly(),

        _ => {
            eprintln!("Invalid parameter: {param:?}");
            return ExitCode::FAILURE;
        }
    }

    ExitCode::SUCCESS
}
