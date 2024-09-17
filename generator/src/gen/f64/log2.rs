use super::{print_f64_const, split_hi_lo};

pub(crate) fn gen_consts() {
    let mut tmp = dev_mpfr::Mpfr::new();
    tmp.set_prec(1024);

    // tmp = log2(e)
    tmp.set_ui(1, dev_mpfr::Rnd::N);
    tmp.exp(None, dev_mpfr::Rnd::N);
    tmp.log2(None, dev_mpfr::Rnd::N);

    let (hi, lo) = split_hi_lo(&mut tmp, 27);

    print_f64_const("LOG2_E_HI", hi);
    print_f64_const("LOG2_E_LO", lo);
}
