use std::io::Write as _;
use std::path::PathBuf;

pub(crate) mod f32;
pub(crate) mod f64;

pub(crate) fn create_prng() -> impl rand::Rng {
    use rand::SeedableRng as _;
    rand_pcg::Pcg64::seed_from_u64(0x985A9231A0046A3D)
}

fn data_path(name: &str) -> PathBuf {
    let dir = env!("CARGO_MANIFEST_DIR");
    let mut path = PathBuf::from(dir);
    path.push("data");
    path.push(name);
    path
}

const HEADER: [u8; 20] = *b"\0\0fpmath\0test\0data\0\0";

fn generate_data<S, T: bincode::Encode>(
    name: &str,
    source: impl FnOnce() -> Vec<S>,
    mut f: impl FnMut(S) -> T,
    pb: indicatif::ProgressBar,
) {
    pb.set_message(String::from(name));

    let pb_sty =
        indicatif::ProgressStyle::with_template("[{elapsed_precise}] {spinner:.cyan/blue} {msg}")
            .unwrap()
            .tick_chars("-\\|/");
    pb.set_style(pb_sty);
    pb.enable_steady_tick(std::time::Duration::from_millis(150));

    let source = source();

    let path = data_path(name);
    let file = std::fs::OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(path)
        .expect("open error");
    let mut file = std::io::BufWriter::new(file);

    bincode::encode_into_std_write(HEADER, &mut file, bincode_config()).expect("encode error");

    bincode::encode_into_std_write(
        u32::try_from(source.len()).unwrap(),
        &mut file,
        bincode_config(),
    )
    .expect("encode error");

    pb.disable_steady_tick();
    pb.set_position(0);
    pb.set_length(source.len().try_into().unwrap());
    let pb_sty = indicatif::ProgressStyle::with_template(
        "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>10}/{len:10} {msg}",
    )
    .unwrap()
    .progress_chars("##-");
    pb.set_style(pb_sty);

    for src_value in source {
        let res_value = f(src_value);
        bincode::encode_into_std_write(res_value, &mut file, bincode_config())
            .expect("encode error");
        pb.inc(1);
    }

    file.flush().expect("write error");

    pb.println(format!("- {name} done"));
    pb.finish_and_clear();
}

#[cfg(test)]
pub(crate) fn consume_data<T: bincode::Decode>(name: &str, mut f: impl FnMut(T)) {
    let path = data_path(name);
    let file = std::fs::OpenOptions::new()
        .read(true)
        .open(path)
        .expect("open error");
    let mut file = std::io::BufReader::new(file);

    let header: [u8; 20] =
        bincode::decode_from_std_read(&mut file, bincode_config()).expect("decode error");
    if header != HEADER {
        panic!("unexpected header");
    }

    let n: u32 = bincode::decode_from_std_read(&mut file, bincode_config()).expect("decode error");
    for _ in 0..n {
        let value =
            bincode::decode_from_std_read(&mut file, bincode_config()).expect("decode error");
        f(value);
    }
}

#[inline]
fn bincode_config() -> impl bincode::config::Config {
    bincode::config::standard()
        .with_little_endian()
        .with_fixed_int_encoding()
}
