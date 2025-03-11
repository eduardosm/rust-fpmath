use crate::RunError;

pub(crate) fn generate(spec: &str) -> Result<String, RunError> {
    let mut splitted = spec.split(' ');
    let cmd = splitted.next().unwrap();
    match cmd {
        "other" => {
            let param = splitted.next().ok_or_else(|| {
                eprintln!("Missing parameter for \"other\" generate command");
                RunError
            })?;
            expect_no_extra(splitted)?;

            let child = std::process::Command::new("cargo")
                .args(["run", "-p", "generator-aux", "--", param])
                .stdin(std::process::Stdio::null())
                .stdout(std::process::Stdio::piped())
                .stderr(std::process::Stdio::inherit())
                .spawn()
                .map_err(|e| {
                    eprintln!("Failed to spawn child: {e}");
                    RunError
                })?;

            let output = child.wait_with_output().map_err(|e| {
                eprintln!("Failed to wait for child: {e}");
                RunError
            })?;
            if !output.status.success() {
                eprintln!("Child returned {}", output.status);
                return Err(RunError);
            }
            assert!(output.stderr.is_empty());

            let stdout = String::from_utf8(output.stdout).map_err(|_| {
                eprintln!("Output is not UTF-8");
                RunError
            })?;

            Ok(stdout)
        }
        _ => {
            eprintln!("Unknown generate command: {cmd:?}");
            Err(RunError)
        }
    }
}

fn expect_no_extra<'a>(mut iter: impl Iterator<Item = &'a str>) -> Result<(), RunError> {
    if let Some(arg) = iter.next() {
        eprintln!("Unexpected extra argument: {arg:?}");
        Err(RunError)
    } else {
        Ok(())
    }
}
