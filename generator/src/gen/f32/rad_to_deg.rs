use super::{print_f32_const, split_hi_lo};

pub(crate) fn gen_consts() {
    let mut tmp = dev_mpfr::Mpfr::new();
    tmp.set_prec(1024);

    // tmp = 180/π
    tmp.const_pi(dev_mpfr::Rnd::N);
    tmp.f64_div(180.0, None, dev_mpfr::Rnd::N);

    let v = tmp.get_f32(dev_mpfr::Rnd::N);
    let (hi, lo) = split_hi_lo(&mut tmp, 12);

    print_f32_const("RAD_TO_DEG", v);
    print_f32_const("RAD_TO_DEG_HI", hi);
    print_f32_const("RAD_TO_DEG_LO", lo);
}
