#![warn(
    rust_2018_idioms,
    trivial_casts,
    trivial_numeric_casts,
    unreachable_pub,
    unused_qualifications
)]

use std::path::{Path, PathBuf};

fn main() {
    let out_dir = PathBuf::from(std::env::var_os("OUT_DIR").unwrap());
    gen_bindings(&out_dir);
}

fn gen_bindings(out_dir: &Path) {
    let header_src = "#include <mpfr.h>\n";
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rustc-link-lib=mpfr");

    let bindings = bindgen::Builder::default()
        .header_contents("mpfr_wrapper.h", header_src)
        .prepend_enum_name(false)
        .default_enum_style(bindgen::EnumVariation::Consts)
        .allowlist_type("mpfr_.*")
        .allowlist_function("mpfr_.*")
        .allowlist_var("MPFR_.*")
        .blocklist_function("mpfr_set_ld|mpfr_get_ld|mpfr_cmp_ld|mpfr_get_ld_2exp")
        .layout_tests(false)
        .rust_target(bindgen::RustTarget::Stable_1_68)
        .formatter(bindgen::Formatter::Rustfmt)
        .generate()
        .unwrap();

    bindings.write_to_file(out_dir.join("mpfr_ffi.rs")).unwrap();
}
