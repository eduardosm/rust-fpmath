use super::{print_f64_const, split_hi_lo};
use crate::sollya;

pub(crate) fn gen_consts() {
    let mut tmp = dev_mpfr::Mpfr::new();
    tmp.set_prec(1024);

    // tmp = cbrt(2)
    tmp.set_ui(2, dev_mpfr::Rnd::N);
    tmp.cbrt(None, dev_mpfr::Rnd::N);

    let (hi, lo) = split_hi_lo(&mut tmp, 27);

    print_f64_const("CBRT_2_HI", hi);
    print_f64_const("CBRT_2_LO", lo);

    // tmp = cbrt(4)
    tmp.set_ui(4, dev_mpfr::Rnd::N);
    tmp.cbrt(None, dev_mpfr::Rnd::N);

    let (hi, lo) = split_hi_lo(&mut tmp, 27);

    print_f64_const("CBRT_4_HI", hi);
    print_f64_const("CBRT_4_LO", lo);
}

pub(crate) fn gen_inv_cbrt_poly() {
    let f = "x^(-1/3)";
    let poly_i = [0, 1, 2, 3, 4, 5];
    let range = (1.0 - 0.001, 2.0 + 0.001);

    sollya::run_and_print_remez_f64(f, range, &poly_i, 0, "K");
}
