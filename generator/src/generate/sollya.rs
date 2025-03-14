use std::fmt::Write as _;
use std::io::{Read as _, Write as _};

use super::FloatKind;

#[derive(Debug)]
enum SollyaError {
    SpawnFailed(std::io::Error),
    StdinWriteFailed(std::io::Error),
    StdoutReadFailed(std::io::Error),
    WaitFailed(std::io::Error),
    ExitError(std::process::ExitStatus),
}

impl std::fmt::Display for SollyaError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::SpawnFailed(e) => write!(f, "failed to spawn sollya: {e}"),
            Self::StdinWriteFailed(e) => write!(f, "failed to write to sollya stdin: {e}"),
            Self::StdoutReadFailed(e) => write!(f, "failed to read from sollya stdout: {e}"),
            Self::WaitFailed(e) => write!(f, "failed to wait for sollya process to finish: {e}"),
            Self::ExitError(status) => write!(f, "sollya process exited with status: {status}"),
        }
    }
}

fn run_sollya(input: &[u8]) -> Result<Vec<u8>, SollyaError> {
    let mut child = std::process::Command::new("sollya")
        .arg("--warnonstderr")
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::inherit())
        .spawn()
        .map_err(SollyaError::SpawnFailed)?;

    let mut child_stdin = child.stdin.take().unwrap();
    let mut child_stdout = child.stdout.take().unwrap();

    let stdout_reader = std::thread::spawn(move || {
        let mut data = Vec::new();
        child_stdout.read_to_end(&mut data)?;
        Ok(data)
    });
    child_stdin
        .write_all(input)
        .map_err(SollyaError::StdinWriteFailed)?;
    drop(child_stdin);

    let stdout_data = stdout_reader
        .join()
        .expect("stdout reader thread panicked")
        .map_err(SollyaError::StdoutReadFailed)?;

    let exit_status = child.wait().map_err(SollyaError::WaitFailed)?;
    if !exit_status.success() {
        return Err(SollyaError::ExitError(exit_status));
    }

    Ok(stdout_data)
}

pub(super) fn run_and_render_remez(
    fkind: FloatKind,
    func: &str,
    range: (f64, f64),
    poly_i: &[i32],
    poly_i_print_off: i32,
    coeff_prefix: &str,
    out: &mut String,
) {
    let code = gen_remez_code(func, range, poly_i);
    let result = run_sollya(&code).unwrap();

    let mut lines = result.split(|&c| c == b'\n');
    let first_line = lines.next().unwrap();
    assert_eq!(first_line, b"The precision has been set to 2048 bits.");

    let err_line = lines.next().unwrap();
    let err = parse_f64(err_line);
    eprintln!("error = {err:e} = 2^({})", err.log2());

    for &i in poly_i.iter() {
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

fn gen_remez_code(func: &str, range: (f64, f64), poly_i: &[i32]) -> Vec<u8> {
    let mut code = Vec::new();
    code.extend(b"prec = 2048;\n");
    code.extend(b"f = ");
    code.extend(func.as_bytes());
    code.extend(b";\n");

    code.extend(b"p = remez(f, [|");
    let mut first = true;
    for &i in poly_i.iter() {
        if !first {
            code.extend(b", ");
        }
        first = false;
        write!(code, "{i}").unwrap();
    }
    writeln!(code, "|], [{}, {}]);", range.0, range.1).unwrap();

    writeln!(
        code,
        "print(dirtyinfnorm(p - f, [{}, {}]));",
        range.0, range.1,
    )
    .unwrap();

    for &i in poly_i.iter() {
        writeln!(code, "print(coeff(p, {i}));").unwrap();
    }
    code.extend(b"quit;\n");

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
