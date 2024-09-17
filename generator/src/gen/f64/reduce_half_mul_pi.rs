use super::{print_f64_const, split_hi_lo};

pub(crate) fn gen_consts() {
    let mut tmp = dev_mpfr::Mpfr::new();
    tmp.set_prec(1024);

    // tmp = π
    tmp.const_pi(dev_mpfr::Rnd::N);

    let v = tmp.get_f64(dev_mpfr::Rnd::N);
    let (hi, lo) = split_hi_lo(&mut tmp, 27);

    print_f64_const("PI", v);
    print_f64_const("PI_HI", hi);
    print_f64_const("PI_LO", lo);
}
