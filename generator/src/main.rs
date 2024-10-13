#![warn(
    rust_2018_idioms,
    trivial_casts,
    trivial_numeric_casts,
    unreachable_pub,
    unused_qualifications
)]
#![forbid(unsafe_code)]

use std::process::ExitCode;

mod gen;
mod julia;
mod sollya;

fn main() -> ExitCode {
    for arg in std::env::args_os().skip(1) {
        println!("{arg:?}:");
        match &*arg.to_string_lossy() {
            // common
            "common::reduce_pi_2_large::FRAC_2_PI_LARGE" => {
                gen::common::reduce_pi_2_large::gen_frac_2_pi_large();
            }
            "common::reduce_pi_2_large::FRAC_PI_2_MEDIUM" => {
                gen::common::reduce_pi_2_large::gen_frac_pi_2_medium();
            }

            // f32
            "f32::consts" => gen::f32::gen_consts(),
            "f32::cbrt::consts" => gen::f32::cbrt::gen_consts(),
            "f32::cbrt::inv_cbrt_poly" => gen::f32::cbrt::gen_inv_cbrt_poly(),
            "f32::div_pi::consts" => gen::f32::div_pi::gen_consts(),
            "f32::exp::consts" => gen::f32::exp::gen_consts(),
            "f32::exp::exp_special_poly" => gen::f32::exp::gen_exp_special_poly(),
            "f32::exp::exp_m1_special_poly" => gen::f32::exp::gen_exp_m1_special_poly(),
            "f32::exp2::consts" => gen::f32::exp2::gen_consts(),
            "f32::exp10::consts" => gen::f32::exp10::gen_consts(),
            "f32::gamma::consts" => gen::f32::gamma::gen_consts(),
            "f32::gamma::lgamma_poly_1" => gen::f32::gamma::gen_lgamma_poly_1(),
            "f32::gamma::lgamma_poly_2" => gen::f32::gamma::gen_lgamma_poly_2(),
            "f32::gamma::special_poly" => gen::f32::gamma::gen_special_poly(),
            "f32::log::consts" => gen::f32::log::gen_consts(),
            "f32::log::log_special_poly" => gen::f32::log::gen_log_special_poly(),
            "f32::log::log_special_poly_ex" => gen::f32::log::gen_log_special_poly_ex(),
            "f32::log2::consts" => gen::f32::log2::gen_consts(),
            "f32::log10::consts" => gen::f32::log10::gen_consts(),
            "f32::rad_to_deg::consts" => gen::f32::rad_to_deg::gen_consts(),
            "f32::reduce_90_deg::consts" => gen::f32::reduce_90_deg::gen_consts(),
            "f32::reduce_half_mul_pi::consts" => gen::f32::reduce_half_mul_pi::gen_consts(),
            "f32::reduce_pi_2::consts" => gen::f32::reduce_pi_2::gen_consts(),
            "f32::sin_cos::consts" => gen::f32::sin_cos::gen_consts(),
            "f32::sin_cos::sin_poly" => gen::f32::sin_cos::gen_sin_poly(),
            "f32::sin_cos::sin_poly_ex" => gen::f32::sin_cos::gen_sin_poly_ex(),
            "f32::sin_cos::cos_poly" => gen::f32::sin_cos::gen_cos_poly(),
            "f32::tan::tan_poly" => gen::f32::tan::gen_tan_poly(),
            "f32::asin_acos::consts" => gen::f32::asin_acos::gen_consts(),
            "f32::asin_acos::asin_poly" => gen::f32::asin_acos::gen_asin_poly(),
            "f32::atan::consts" => gen::f32::atan::gen_consts(),
            "f32::atan::atan_poly" => gen::f32::atan::gen_atan_poly(),

            // f64
            "f64::consts" => gen::f64::gen_consts(),
            "f64::cbrt::consts" => gen::f64::cbrt::gen_consts(),
            "f64::cbrt::inv_cbrt_poly" => gen::f64::cbrt::gen_inv_cbrt_poly(),
            "f64::div_pi::consts" => gen::f64::div_pi::gen_consts(),
            "f64::exp::consts" => gen::f64::exp::gen_consts(),
            "f64::exp::exp_special_poly" => gen::f64::exp::gen_exp_special_poly(),
            "f64::exp::exp_m1_special_poly" => gen::f64::exp::gen_exp_m1_special_poly(),
            "f64::exp2::consts" => gen::f64::exp2::gen_consts(),
            "f64::exp10::consts" => gen::f64::exp10::gen_consts(),
            "f64::gamma::consts" => gen::f64::gamma::gen_consts(),
            "f64::gamma::lgamma_poly_1" => gen::f64::gamma::gen_lgamma_poly_1(),
            "f64::gamma::lgamma_poly_2" => gen::f64::gamma::gen_lgamma_poly_2(),
            "f64::gamma::special_poly" => gen::f64::gamma::gen_special_poly(),
            "f64::log::consts" => gen::f64::log::gen_consts(),
            "f64::log::log_special_poly" => gen::f64::log::gen_log_special_poly(),
            "f64::log::log_special_poly_ex" => gen::f64::log::gen_log_special_poly_ex(),
            "f64::log2::consts" => gen::f64::log2::gen_consts(),
            "f64::log10::consts" => gen::f64::log10::gen_consts(),
            "f64::rad_to_deg::consts" => gen::f64::rad_to_deg::gen_consts(),
            "f64::reduce_90_deg::consts" => gen::f64::reduce_90_deg::gen_consts(),
            "f64::reduce_half_mul_pi::consts" => gen::f64::reduce_half_mul_pi::gen_consts(),
            "f64::reduce_pi_2::consts" => gen::f64::reduce_pi_2::gen_consts(),
            "f64::sin_cos::consts" => gen::f64::sin_cos::gen_consts(),
            "f64::sin_cos::sin_poly" => gen::f64::sin_cos::gen_sin_poly(),
            "f64::sin_cos::sin_poly_ex" => gen::f64::sin_cos::gen_sin_poly_ex(),
            "f64::sin_cos::cos_poly" => gen::f64::sin_cos::gen_cos_poly(),
            "f64::tan::tan_poly" => gen::f64::tan::gen_tan_poly(),
            "f64::asin_acos::consts" => gen::f64::asin_acos::gen_consts(),
            "f64::asin_acos::asin_poly" => gen::f64::asin_acos::gen_asin_poly(),
            "f64::atan::consts" => gen::f64::atan::gen_consts(),
            "f64::atan::atan_poly" => gen::f64::atan::gen_atan_poly(),

            arg => {
                eprintln!("invalid argument {arg:?}");
                return ExitCode::FAILURE;
            }
        }
    }

    ExitCode::SUCCESS
}
