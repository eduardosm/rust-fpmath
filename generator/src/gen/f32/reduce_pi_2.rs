use super::{print_f32_const, split_hi_lo};

pub(crate) fn gen_consts() {
    let mut tmp = dev_mpfr::Mpfr::new();
    tmp.set_prec(1024);

    // tmp = π/2
    tmp.const_pi(dev_mpfr::Rnd::N);
    tmp.div_f64(None, 2.0, dev_mpfr::Rnd::N);

    let v = tmp.get_f32(dev_mpfr::Rnd::N);

    print_f32_const("FRAC_PI_2", v);

    // tmp = π/4
    tmp.const_pi(dev_mpfr::Rnd::N);
    tmp.div_f64(None, 4.0, dev_mpfr::Rnd::N);

    let v = tmp.get_f32(dev_mpfr::Rnd::N);

    print_f32_const("FRAC_PI_4", v);

    // tmp = 2/π
    tmp.const_pi(dev_mpfr::Rnd::N);
    tmp.f64_div(2.0, None, dev_mpfr::Rnd::N);

    let v = tmp.get_f32(dev_mpfr::Rnd::N);

    print_f32_const("FRAC_2_PI", v);

    // tmp = π/2
    tmp.const_pi(dev_mpfr::Rnd::N);
    tmp.mul_f64(None, 0.5, dev_mpfr::Rnd::N);

    let (hi, hiex) = split_hi_lo(&mut tmp, 9);
    let (mi, miex) = split_hi_lo(&mut tmp, 9);
    let (lo, loex) = split_hi_lo(&mut tmp, 9);

    print_f32_const("FRAC_PI_2_HI", hi);
    print_f32_const("FRAC_PI_2_HIEX", hiex);
    print_f32_const("FRAC_PI_2_MI", mi);
    print_f32_const("FRAC_PI_2_MIEX", miex);
    print_f32_const("FRAC_PI_2_LO", lo);
    print_f32_const("FRAC_PI_2_LOEX", loex);
}
