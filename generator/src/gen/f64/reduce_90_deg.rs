use super::{print_f64_const, split_hi_lo};

pub(crate) fn gen_consts() {
    let mut tmp = dev_mpfr::Mpfr::new();
    tmp.set_prec(1024);

    // tmp = π/180
    tmp.const_pi(dev_mpfr::Rnd::N);
    tmp.div_f64(None, 180.0, dev_mpfr::Rnd::N);

    let v = tmp.get_f64(dev_mpfr::Rnd::N);
    let (hi, lo) = split_hi_lo(&mut tmp, 27);

    print_f64_const("DEG_TO_RAD", v);
    print_f64_const("DEG_TO_RAD_HI", hi);
    print_f64_const("DEG_TO_RAD_LO", lo);
}
