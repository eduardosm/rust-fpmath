#![warn(
    rust_2018_idioms,
    trivial_casts,
    trivial_numeric_casts,
    unreachable_pub,
    unused_qualifications
)]
#![forbid(unsafe_code)]
#![allow(clippy::too_many_arguments)]

use std::ffi::OsStr;
use std::path::Path;
use std::process::ExitCode;

mod generate;

struct RunError;

fn main() -> ExitCode {
    let mut args = std::env::args_os();
    let arg0 = args.next().unwrap();

    if args.len() == 0 {
        eprintln!("Usage: {} <paths...>", arg0.to_string_lossy());
        return ExitCode::FAILURE;
    }

    for arg in args {
        if let Err(RunError) = proc_path(Path::new(&arg)) {
            return ExitCode::FAILURE;
        }
    }

    ExitCode::SUCCESS
}

fn proc_path(path: &Path) -> Result<(), RunError> {
    let meta = path.metadata().map_err(|e| {
        eprintln!("Failed to get metadata of {path:?}: {e}");
        RunError
    })?;

    if meta.is_dir() {
        proc_dir(path)
    } else {
        proc_file(path)
    }
}

fn proc_dir(path: &Path) -> Result<(), RunError> {
    let dir = path.read_dir().map_err(|e| {
        eprintln!("Failed to read directory {path:?}: {e}");
        RunError
    })?;
    for entry in dir {
        let entry = entry.map_err(|e| {
            eprintln!("Failed to read directory {path:?}: {e}");
            RunError
        })?;
        proc_path(&entry.path())?;
    }
    Ok(())
}

fn proc_file(path: &Path) -> Result<(), RunError> {
    if path.extension() != Some(OsStr::new("rs")) {
        return Ok(());
    }

    eprintln!("Processing file {path:?}");

    let file_data = std::fs::read_to_string(path).map_err(|e| {
        eprintln!("Failed to read file {path:?}: {e}");
        RunError
    })?;
    let mut result = String::new();

    let mut num_generates = 0usize;
    let mut lines = file_data.lines();
    while let Some(line) = lines.next() {
        result.push_str(line);
        result.push('\n');
        if let Some((ident, txt)) = check_comment(line, check_generate) {
            let end_line = lines
                .find_map(|line| {
                    let line_ident_len = line.find(|c| c != ' ').unwrap_or(0);
                    if line_ident_len < ident.len() || line.trim_matches(' ').is_empty() {
                        Some(Ok(line))
                    } else if check_comment(line, check_generate).is_some() {
                        eprintln!("Another \"GENERATE\" found before expected");
                        Some(Err(RunError))
                    } else {
                        None
                    }
                })
                .ok_or_else(|| {
                    eprintln!("End line not found for \"GENERATE\" block");
                    RunError
                })??;
            eprintln!("GENERATE: {}", txt.trim());
            let generated = generate::generate(txt)?;
            for generated_line in generated.lines() {
                result.push_str(ident);
                result.push_str(generated_line);
                result.push('\n');
            }
            result.push_str(end_line);
            result.push('\n');
            num_generates += 1;
        }
    }

    std::fs::write(path, result).map_err(|e| {
        eprintln!("Failed to write file {path:?}: {e}");
        RunError
    })?;

    eprintln!("Processed {num_generates} \"GENERATE\" blocks");

    Ok(())
}

fn check_comment<'a, R: 'a>(
    line: &'a str,
    f: impl FnOnce(&'a str) -> Option<R>,
) -> Option<(&'a str, R)> {
    let no_sp_i = line.find(|c| c != ' ')?;
    let (ident, txt) = line.split_at(no_sp_i);
    let txt = txt.strip_prefix("//")?.trim_start_matches(' ');
    f(txt).map(|r| (ident, r))
}

fn check_generate(txt: &str) -> Option<&str> {
    txt.strip_prefix("GENERATE:")
}
