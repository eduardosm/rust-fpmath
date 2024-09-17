use super::{print_f64_const, split_hi_lo};
use crate::sollya;

pub(crate) fn gen_consts() {
    let mut tmp = dev_mpfr::Mpfr::new();
    tmp.set_prec(1024);

    // tmp = sqrt(2)
    tmp.set_ui(2, dev_mpfr::Rnd::N);
    tmp.sqrt(None, dev_mpfr::Rnd::N);

    let v = tmp.get_f64(dev_mpfr::Rnd::N);

    print_f64_const("SQRT_2", v);

    // tmp = ln(2)
    tmp.set_ui(2, dev_mpfr::Rnd::N);
    tmp.log(None, dev_mpfr::Rnd::N);

    let (hi, lo) = split_hi_lo(&mut tmp, 27);

    print_f64_const("LN_2_HI", hi);
    print_f64_const("LN_2_LO", lo);

    // tmp = 2 / 3
    tmp.set_ui(2, dev_mpfr::Rnd::N);
    tmp.div_f64(None, 3.0, dev_mpfr::Rnd::N);

    let (hi, lo) = split_hi_lo(&mut tmp, 0);

    print_f64_const("FRAC_2_3_HI", hi);
    print_f64_const("FRAC_2_3_LO", lo);

    // tmp = 0.4
    tmp.set_ui(4, dev_mpfr::Rnd::N);
    tmp.div_f64(None, 10.0, dev_mpfr::Rnd::N);

    let (hi, lo) = split_hi_lo(&mut tmp, 0);

    print_f64_const("FRAC_4_10_HI", hi);
    print_f64_const("FRAC_4_10_LO", lo);
}

pub(crate) fn gen_log_special_poly() {
    let f = "(log1p(x) - log(1 - x) - 2 * x) / x";
    let poly_i = [2, 4, 6, 8, 10, 12, 14, 99];
    let range0 = 0.1716;
    let range = (-range0, range0);

    sollya::run_and_print_remez_f64(f, range, &poly_i, 0, "K");
}

pub(crate) fn gen_log_special_poly_ex() {
    let f = "(log1p(x) - log(1 - x) - 2 * x - (2/3) * x^3 - 0.4 * x^5) / x";
    let poly_i = [6, 8, 10, 12, 14, 16];
    let range = (-0.1, 0.1716);

    sollya::run_and_print_remez_f64(f, range, &poly_i, 0, "K");
}
