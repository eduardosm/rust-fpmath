use super::{print_f64_const, split_hi_lo};

pub(crate) fn gen_consts() {
    let mut tmp = dev_mpfr::Mpfr::new();
    tmp.set_prec(1024);

    // tmp = ln(2)
    tmp.set_ui(2, dev_mpfr::Rnd::N);
    tmp.log(None, dev_mpfr::Rnd::N);

    let v = tmp.get_f64(dev_mpfr::Rnd::N);
    let (hi, lo) = split_hi_lo(&mut tmp, 27);

    print_f64_const("LN_2", v);
    print_f64_const("LN_2_HI", hi);
    print_f64_const("LN_2_LO", lo);
}
