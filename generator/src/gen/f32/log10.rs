use super::{print_f32_const, split_hi_lo};

pub(crate) fn gen_consts() {
    let mut tmp = dev_mpfr::Mpfr::new();
    tmp.set_prec(1024);

    // tmp = log10(e)
    tmp.set_ui(1, dev_mpfr::Rnd::N);
    tmp.exp(None, dev_mpfr::Rnd::N);
    tmp.log10(None, dev_mpfr::Rnd::N);

    let (hi, lo) = split_hi_lo(&mut tmp, 12);

    print_f32_const("LOG10_E_HI", hi);
    print_f32_const("LOG10_E_LO", lo);

    // tmp = log10(2)
    tmp.set_ui(2, dev_mpfr::Rnd::N);
    tmp.log10(None, dev_mpfr::Rnd::N);

    let (hi, lo) = split_hi_lo(&mut tmp, 12);

    print_f32_const("LOG10_2_HI", hi);
    print_f32_const("LOG10_2_LO", lo);
}
