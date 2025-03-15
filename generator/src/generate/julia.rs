use std::fmt::Write as _;
use std::io::Read as _;

use super::FloatKind;

#[derive(Debug)]
enum JuliaError {
    SpawnFailed(std::io::Error),
    StdoutReadFailed(std::io::Error),
    WaitFailed(std::io::Error),
    ExitError(std::process::ExitStatus),
}

impl std::fmt::Display for JuliaError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::SpawnFailed(e) => write!(f, "failed to spawn julia: {e}"),
            Self::StdoutReadFailed(e) => write!(f, "failed to read from julia stdout: {e}"),
            Self::WaitFailed(e) => write!(f, "failed to wait for julia process to finish: {e}"),
            Self::ExitError(status) => write!(f, "julia process exited with status: {status}"),
        }
    }
}

fn run_julia(input: &str) -> Result<Vec<u8>, JuliaError> {
    let mut child = std::process::Command::new("julia")
        .arg("-e")
        .arg(input)
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::inherit())
        .spawn()
        .map_err(JuliaError::SpawnFailed)?;

    let mut child_stdout = child.stdout.take().unwrap();

    let mut stdout_data = Vec::new();
    child_stdout
        .read_to_end(&mut stdout_data)
        .map_err(JuliaError::StdoutReadFailed)?;

    let exit_status = child.wait().map_err(JuliaError::WaitFailed)?;
    if !exit_status.success() {
        return Err(JuliaError::ExitError(exit_status));
    }

    Ok(stdout_data)
}

pub(super) fn run_and_render_remez(
    fkind: FloatKind,
    func: &str,
    wfunc: &str,
    range: (f64, f64),
    poly_deg: i32,
    poly_i_print_off: i32,
    coeff_prefix: &str,
    out: &mut String,
) {
    let code = gen_remez_code(func, wfunc, range, poly_deg);
    let result = run_julia(&code).unwrap();

    let mut lines = result.split(|&c| c == b'\n');

    let err_line = lines.next().unwrap();
    let err = parse_f64(err_line);
    eprintln!("error = {err:e} = 2^({})", err.log2());

    for i in 0..=poly_deg {
        let coeff_line = lines.next().unwrap();
        match fkind {
            FloatKind::F32 => {
                let coeff_value = parse_f32(coeff_line);
                writeln!(
                    out,
                    "const {coeff_prefix}{}: u32 = 0x{:08X}; // {coeff_value:e}",
                    i + poly_i_print_off,
                    coeff_value.to_bits(),
                )
                .unwrap();
            }
            FloatKind::F64 => {
                let coeff_value = parse_f64(coeff_line);
                writeln!(
                    out,
                    "const {coeff_prefix}{}: u64 = 0x{:016X}; // {coeff_value:e}",
                    i + poly_i_print_off,
                    coeff_value.to_bits(),
                )
                .unwrap();
            }
        }
    }
}

fn gen_remez_code(func: &str, wfunc: &str, range: (f64, f64), poly_deg: i32) -> String {
    let mut code = String::new();

    code.push_str("import Remez;\n");
    code.push_str("import SpecialFunctions;\n");
    code.push_str("Remez.setprecision(BigFloat, 1024);\n");
    code.push_str("f = (x) -> ");
    code.push_str(func);
    code.push_str(";\n");
    code.push_str("w = (x, fx) -> ");
    code.push_str(wfunc);
    code.push_str(";\n");

    code.push_str("N, D, E, X = Remez.ratfn_minimax(f, ");
    write!(code, "({}, {}), ", range.0, range.1).unwrap();
    write!(code, "{poly_deg},").unwrap();
    code.push_str(" 0, w);\n");

    code.push_str("println(E);\n");
    for i in 1..=(poly_deg + 1) {
        writeln!(code, "println(N[{i}]);",).unwrap();
    }

    code
}

fn parse_f32(s: &[u8]) -> f32 {
    let s = std::str::from_utf8(s).unwrap();
    s.parse().unwrap()
}

fn parse_f64(s: &[u8]) -> f64 {
    let s = std::str::from_utf8(s).unwrap();
    s.parse().unwrap()
}
