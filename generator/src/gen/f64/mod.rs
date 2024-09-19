pub(crate) mod asin_acos;
pub(crate) mod atan;
pub(crate) mod cbrt;
pub(crate) mod div_pi;
pub(crate) mod exp;
pub(crate) mod exp10;
pub(crate) mod exp2;
pub(crate) mod log;
pub(crate) mod log10;
pub(crate) mod log2;
pub(crate) mod rad_to_deg;
pub(crate) mod reduce_90_deg;
pub(crate) mod reduce_half_mul_pi;
pub(crate) mod reduce_pi_2;
pub(crate) mod sin_cos;
pub(crate) mod tan;

fn split_hi_lo(tmp: &mut dev_mpfr::Mpfr, n_zeros_in_hi: u8) -> (f64, f64) {
    let hi = tmp.get_f64(dev_mpfr::Rnd::Z);
    let hi = f64::from_bits(hi.to_bits() & (u64::MAX << n_zeros_in_hi));
    tmp.sub_f64(None, hi, dev_mpfr::Rnd::N);
    let lo = tmp.get_f64(dev_mpfr::Rnd::N);
    (hi, lo)
}

fn print_f64_const<N: std::fmt::Display>(name: N, value: f64) {
    println!(
        "const {name}: u64 = 0x{:016X}; // {value:e}",
        value.to_bits(),
    );
}

pub(crate) fn gen_consts() {
    let mut tmp = dev_mpfr::Mpfr::new();
    tmp.set_prec(1024);

    // tmp = π
    tmp.const_pi(dev_mpfr::Rnd::N);

    let v = tmp.get_f64(dev_mpfr::Rnd::N);
    print_f64_const("PI", v);

    // tmp = π/2
    tmp.const_pi(dev_mpfr::Rnd::N);
    tmp.div_f64(None, 2.0, dev_mpfr::Rnd::N);

    let v = tmp.get_f64(dev_mpfr::Rnd::N);
    print_f64_const("FRAC_PI_2", v);

    // tmp = π/4
    tmp.const_pi(dev_mpfr::Rnd::N);
    tmp.div_f64(None, 4.0, dev_mpfr::Rnd::N);

    let v = tmp.get_f64(dev_mpfr::Rnd::N);
    print_f64_const("FRAC_PI_4", v);

    // tmp = 1/π
    tmp.const_pi(dev_mpfr::Rnd::N);
    tmp.f64_div(1.0, None, dev_mpfr::Rnd::N);

    let v = tmp.get_f64(dev_mpfr::Rnd::N);
    print_f64_const("FRAC_1_PI", v);

    // tmp = 2/π
    tmp.const_pi(dev_mpfr::Rnd::N);
    tmp.f64_div(2.0, None, dev_mpfr::Rnd::N);

    let v = tmp.get_f64(dev_mpfr::Rnd::N);
    print_f64_const("FRAC_2_PI", v);
}
