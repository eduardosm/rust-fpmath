#![warn(
    rust_2018_idioms,
    trivial_casts,
    trivial_numeric_casts,
    unreachable_pub,
    unused_qualifications
)]
#![allow(clippy::type_complexity)]
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
            if arg == "--all" {
                tasks.extend(data::GEN_FUNCTIONS.iter().map(|&(_, f)| f));
            } else {
                match data::GEN_FUNCTIONS.iter().find(|(name, _)| arg == name) {
                    Some(&(_, f)) => tasks.push(f),
                    None => {
                        eprintln!("invalid argument {arg:?}");
                        return ExitCode::FAILURE;
                    }
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
