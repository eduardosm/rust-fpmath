use super::{print_f32_const, split_hi_lo};

pub(crate) fn gen_consts() {
    let mut tmp = dev_mpfr::Mpfr::new();
    tmp.set_prec(1024);

    // tmp = log2(10)
    tmp.set_ui(10, dev_mpfr::Rnd::N);
    tmp.log2(None, dev_mpfr::Rnd::N);

    let v = tmp.get_f32(dev_mpfr::Rnd::N);

    print_f32_const("LOG2_10", v);

    // tmp = log10(2)
    tmp.set_ui(2, dev_mpfr::Rnd::N);
    tmp.log10(None, dev_mpfr::Rnd::N);

    let (hi, lo) = split_hi_lo(&mut tmp, 12);

    print_f32_const("LOG10_2_HI", hi);
    print_f32_const("LOG10_2_LO", lo);

    // tmp = ln(10)
    tmp.set_ui(10, dev_mpfr::Rnd::N);
    tmp.log(None, dev_mpfr::Rnd::N);

    let v = tmp.get_f32(dev_mpfr::Rnd::N);
    let (hi, lo) = split_hi_lo(&mut tmp, 12);

    print_f32_const("LN_10", v);
    print_f32_const("LN_10_HI", hi);
    print_f32_const("LN_10_LO", lo);
}
