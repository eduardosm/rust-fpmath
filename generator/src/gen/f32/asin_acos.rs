use super::{print_f32_const, split_hi_lo};
use crate::sollya;

pub(crate) fn gen_consts() {
    let mut tmp = dev_mpfr::Mpfr::new();
    tmp.set_prec(1024);

    // tmp = π/2
    tmp.const_pi(dev_mpfr::Rnd::N);
    tmp.div_f64(None, 2.0, dev_mpfr::Rnd::N);

    let (hi, lo) = split_hi_lo(&mut tmp, 0);

    print_f32_const("FRAC_PI_2_HI", hi);
    print_f32_const("FRAC_PI_2_LO", lo);
}

pub(crate) fn gen_asin_poly() {
    let f = "asin(x) - x";
    let poly_i = [3, 5, 7, 9, 11];
    let range = (-0.001, 0.501);

    sollya::run_and_print_remez_f32(f, range, &poly_i, -3, "K");
}
