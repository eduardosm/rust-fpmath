use std::fmt::Write as _;
use std::io::Read as _;

#[derive(Debug)]
pub(crate) enum JuliaError {
    SpawnFailed(std::io::Error),
    StdoutReadFailed(std::io::Error),
    WaitFailed(std::io::Error),
    ExitError(std::process::ExitStatus),
}

impl std::fmt::Display for JuliaError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::SpawnFailed(e) => write!(f, "failed to spawn sollya: {e}"),
            Self::StdoutReadFailed(e) => write!(f, "failed to read from sollya stdout: {e}"),
            Self::WaitFailed(e) => write!(f, "failed to wait for sollya process to finish: {e}"),
            Self::ExitError(status) => write!(f, "sollya process exited with status: {status}"),
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

pub(crate) fn run_and_print_remez_f32(
    f: &str,
    range: (f64, f64),
    poly_deg: i32,
    poly_i_print_off: i32,
    coeff_prefix: &str,
    exclude_coeffs: &[i32],
) {
    let code = gen_remez_code(f, range, poly_deg, "Float32", "UInt32", 8);
    let result = run_julia(&code).unwrap();

    let mut lines = result.split(|&c| c == b'\n');

    let err_line = lines.next().unwrap();
    let err = parse_f64_hex(err_line);
    eprintln!("error = {err:e} = 2^({})", err.log2());

    for i in 0..=poly_deg {
        let coeff_line = lines.next().unwrap();
        let coeff_value = parse_f32_hex(coeff_line);
        if !exclude_coeffs.contains(&i) {
            println!(
                "const {coeff_prefix}{}: u32 = 0x{:08X}; // {coeff_value:e}",
                i + poly_i_print_off,
                coeff_value.to_bits(),
            );
        }
    }
}

pub(crate) fn run_and_print_remez_f64(
    f: &str,
    range: (f64, f64),
    poly_deg: i32,
    poly_i_print_off: i32,
    coeff_prefix: &str,
    exclude_coeffs: &[i32],
) {
    let code = gen_remez_code(f, range, poly_deg, "Float64", "UInt64", 16);
    let result = run_julia(&code).unwrap();

    let mut lines = result.split(|&c| c == b'\n');

    let err_line = lines.next().unwrap();
    let err = parse_f64_hex(err_line);
    eprintln!("error = {err:e} = 2^({})", err.log2());

    for i in 0..=poly_deg {
        let coeff_line = lines.next().unwrap();
        let coeff_value = parse_f64_hex(coeff_line);
        if !exclude_coeffs.contains(&i) {
            println!(
                "const {coeff_prefix}{}: u64 = 0x{:016X}; // {coeff_value:e}",
                i + poly_i_print_off,
                coeff_value.to_bits(),
            );
        }
    }
}

fn gen_remez_code(
    f: &str,
    range: (f64, f64),
    poly_deg: i32,
    float_type: &str,
    float_bits_type: &str,
    float_nibbles: u8,
) -> String {
    let mut code = String::new();

    code.push_str("import Remez;\n");
    code.push_str("import Printf;\n");
    code.push_str("import SpecialFunctions;\n");
    code.push_str("Remez.setprecision(BigFloat, 512);\n");
    code.push_str("f = (x) -> ");
    code.push_str(f);
    code.push_str(";\n");

    code.push_str("N, D, E, X = Remez.ratfn_minimax(f, ");
    write!(code, "({}, {}), ", range.0, range.1).unwrap();
    write!(code, "{poly_deg},").unwrap();
    code.push_str(" 0);\n");

    code.push_str("Printf.@printf \"0x%016x\\n\" reinterpret(UInt64, Float64(E))\n");
    for i in 1..=(poly_deg + 1) {
        writeln!(
            code,
            "Printf.@printf \"0x%0{float_nibbles}x\\n\" reinterpret({float_bits_type}, {float_type}(N[{i}]));",
        )
        .unwrap();
    }

    code
}

fn parse_f32_hex(s: &[u8]) -> f32 {
    assert!(s.len() == 10 && s.starts_with(b"0x"));
    let mut x: u32 = 0;
    for &c in s[2..].iter() {
        x <<= 4;
        match c {
            b'0'..=b'9' => x += u32::from(c - b'0'),
            b'A'..=b'F' => x += u32::from(c - b'A' + 10),
            b'a'..=b'f' => x += u32::from(c - b'a' + 10),
            _ => panic!("invalid hex digit"),
        }
    }
    f32::from_bits(x)
}

fn parse_f64_hex(s: &[u8]) -> f64 {
    assert!(s.len() == 18 && s.starts_with(b"0x"));
    let mut x: u64 = 0;
    for &c in s[2..].iter() {
        x <<= 4;
        match c {
            b'0'..=b'9' => x += u64::from(c - b'0'),
            b'A'..=b'F' => x += u64::from(c - b'A' + 10),
            b'a'..=b'f' => x += u64::from(c - b'a' + 10),
            _ => panic!("invalid hex digit"),
        }
    }
    f64::from_bits(x)
}
