#![warn(
    rust_2018_idioms,
    trivial_casts,
    trivial_numeric_casts,
    unreachable_pub,
    unused_qualifications
)]
#![forbid(unsafe_code)]

use std::collections::{BTreeSet, HashMap};
use std::process::ExitCode;

mod data;
#[cfg(test)]
mod tests;

fn main() -> ExitCode {
    let mut args = std::env::args_os();
    args.next().unwrap();

    if args.next().map_or(false, |arg| arg == "--generate") {
        let mp = indicatif::MultiProgress::new();

        let mut args_dedup = BTreeSet::new();
        for arg in args {
            match arg.into_string() {
                Ok(arg) => {
                    args_dedup.insert(arg);
                }
                Err(arg) => {
                    eprintln!("invalid argument {arg:?}");
                    return ExitCode::FAILURE;
                }
            }
        }

        let mut tasks = Vec::<fn(indicatif::ProgressBar)>::new();
        for arg in args_dedup.iter() {
            match arg.as_str() {
                "f32::acosh" => tasks.push(data::f32::acosh::gen_data),
                "f32::asin_acos" => tasks.push(data::f32::asin_acos::gen_data),
                "f32::asind_acosd" => tasks.push(data::f32::asin_acos::gen_data_d),
                "f32::asinh" => tasks.push(data::f32::asinh::gen_data),
                "f32::asinpi_acospi" => tasks.push(data::f32::asin_acos::gen_data_pi),
                "f32::atan" => tasks.push(data::f32::atan::gen_data),
                "f32::atan2" => tasks.push(data::f32::atan2::gen_data),
                "f32::atan2d" => tasks.push(data::f32::atan2::gen_data_d),
                "f32::atan2pi" => tasks.push(data::f32::atan2::gen_data_pi),
                "f32::atand" => tasks.push(data::f32::atan::gen_data_d),
                "f32::atanh" => tasks.push(data::f32::atanh::gen_data),
                "f32::atanpi" => tasks.push(data::f32::atan::gen_data_pi),
                "f32::cbrt" => tasks.push(data::f32::cbrt::gen_data),
                "f32::exp" => tasks.push(data::f32::exp::gen_data),
                "f32::exp10" => tasks.push(data::f32::exp::gen_data_10),
                "f32::exp2" => tasks.push(data::f32::exp::gen_data_2),
                "f32::hypot" => tasks.push(data::f32::hypot::gen_data),
                "f32::log" => tasks.push(data::f32::log::gen_data),
                "f32::log10" => tasks.push(data::f32::log::gen_data_10),
                "f32::log2" => tasks.push(data::f32::log::gen_data_2),
                "f32::log_1p" => tasks.push(data::f32::log_1p::gen_data),
                "f32::pow" => tasks.push(data::f32::pow::gen_data),
                "f32::powi" => tasks.push(data::f32::powi::gen_data),
                "f32::sin_cos" => tasks.push(data::f32::sin_cos::gen_data),
                "f32::sind_cosd" => tasks.push(data::f32::sind_cosd::gen_data),
                "f32::sinh_cosh" => tasks.push(data::f32::sinh_cosh::gen_data),
                "f32::sinpi_cospi" => tasks.push(data::f32::sinpi_cospi::gen_data),
                "f32::sqrt" => tasks.push(data::f32::sqrt::gen_data),
                "f32::tan" => tasks.push(data::f32::tan::gen_data),
                "f32::tand" => tasks.push(data::f32::tand::gen_data),
                "f32::tanh" => tasks.push(data::f32::tanh::gen_data),
                "f32::tanpi" => tasks.push(data::f32::tanpi::gen_data),
                "f64::acosh" => tasks.push(data::f64::acosh::gen_data),
                "f64::asin_acos" => tasks.push(data::f64::asin_acos::gen_data),
                "f64::asind_acosd" => tasks.push(data::f64::asin_acos::gen_data_d),
                "f64::asinh" => tasks.push(data::f64::asinh::gen_data),
                "f64::asinpi_acospi" => tasks.push(data::f64::asin_acos::gen_data_pi),
                "f64::atan" => tasks.push(data::f64::atan::gen_data),
                "f64::atan2" => tasks.push(data::f64::atan2::gen_data),
                "f64::atan2d" => tasks.push(data::f64::atan2::gen_data_d),
                "f64::atan2pi" => tasks.push(data::f64::atan2::gen_data_pi),
                "f64::atand" => tasks.push(data::f64::atan::gen_data_d),
                "f64::atanh" => tasks.push(data::f64::atanh::gen_data),
                "f64::atanpi" => tasks.push(data::f64::atan::gen_data_pi),
                "f64::cbrt" => tasks.push(data::f64::cbrt::gen_data),
                "f64::exp" => tasks.push(data::f64::exp::gen_data),
                "f64::exp10" => tasks.push(data::f64::exp::gen_data_10),
                "f64::exp2" => tasks.push(data::f64::exp::gen_data_2),
                "f64::hypot" => tasks.push(data::f64::hypot::gen_data),
                "f64::log" => tasks.push(data::f64::log::gen_data),
                "f64::log10" => tasks.push(data::f64::log::gen_data_10),
                "f64::log2" => tasks.push(data::f64::log::gen_data_2),
                "f64::log_1p" => tasks.push(data::f64::log_1p::gen_data),
                "f64::pow" => tasks.push(data::f64::pow::gen_data),
                "f64::powi" => tasks.push(data::f64::powi::gen_data),
                "f64::sin_cos" => tasks.push(data::f64::sin_cos::gen_data),
                "f64::sind_cosd" => tasks.push(data::f64::sind_cosd::gen_data),
                "f64::sinh_cosh" => tasks.push(data::f64::sinh_cosh::gen_data),
                "f64::sinpi_cospi" => tasks.push(data::f64::sinpi_cospi::gen_data),
                "f64::sqrt" => tasks.push(data::f64::sqrt::gen_data),
                "f64::tan" => tasks.push(data::f64::tan::gen_data),
                "f64::tand" => tasks.push(data::f64::tand::gen_data),
                "f64::tanh" => tasks.push(data::f64::tanh::gen_data),
                "f64::tanpi" => tasks.push(data::f64::tanpi::gen_data),
                _ => {
                    eprintln!("invalid argument {arg:?}");
                    return ExitCode::FAILURE;
                }
            }
        }

        let mut task_iter = tasks.iter();

        let max_threads = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(4);
        let mut running_tasks = HashMap::new();
        let (sender, receiver) = std::sync::mpsc::sync_channel(max_threads);

        loop {
            while running_tasks.len() < max_threads {
                if let Some(&task) = task_iter.next() {
                    let pb = mp.add(indicatif::ProgressBar::new(0));
                    let sender_clone = sender.clone();
                    let join_handle = std::thread::spawn(move || {
                        task(pb);
                        sender_clone.send(std::thread::current().id()).unwrap();
                    });
                    let thread_id = join_handle.thread().id();
                    running_tasks.insert(thread_id, join_handle);
                } else {
                    break;
                }
            }

            if running_tasks.is_empty() {
                break;
            }

            let finished_id = receiver.recv().unwrap();
            running_tasks.remove(&finished_id).unwrap().join().unwrap();
            loop {
                match receiver.try_recv() {
                    Ok(finished_id) => running_tasks.remove(&finished_id).unwrap().join().unwrap(),
                    Err(std::sync::mpsc::TryRecvError::Empty) => break,
                    Err(e) => panic!("{e:?}"),
                }
            }
        }

        ExitCode::SUCCESS
    } else {
        eprintln!("invalid arguments");
        ExitCode::FAILURE
    }
}
