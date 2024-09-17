use super::{print_f32_const, split_hi_lo};
use crate::sollya;

pub(crate) fn gen_consts() {
    let mut tmp = dev_mpfr::Mpfr::new();
    tmp.set_prec(1024);

    // tmp = π
    tmp.const_pi(dev_mpfr::Rnd::N);

    let v = tmp.get_f32(dev_mpfr::Rnd::N);

    print_f32_const("PI", v);

    // tmp = π/2
    tmp.const_pi(dev_mpfr::Rnd::N);
    tmp.div_f64(None, 2.0, dev_mpfr::Rnd::N);

    let v = tmp.get_f32(dev_mpfr::Rnd::N);
    let (hi, lo) = split_hi_lo(&mut tmp, 0);

    print_f32_const("FRAC_PI_2", v);
    print_f32_const("FRAC_PI_2_HI", hi);
    print_f32_const("FRAC_PI_2_LO", lo);

    // tmp = π/4
    tmp.const_pi(dev_mpfr::Rnd::N);
    tmp.div_f64(None, 4.0, dev_mpfr::Rnd::N);

    let v = tmp.get_f32(dev_mpfr::Rnd::N);

    print_f32_const("FRAC_PI_4", v);

    // tmp = 3π/4
    tmp.const_pi(dev_mpfr::Rnd::N);
    tmp.mul_f64(None, 0.75, dev_mpfr::Rnd::N);

    let v = tmp.get_f32(dev_mpfr::Rnd::N);

    print_f32_const("FRAC_3PI_4", v);
}

pub(crate) fn gen_atan_poly() {
    let f = "atan(x) - x";
    let poly_i = [3, 5, 7, 9, 11, 13, 15, 17, 19];
    let range = (-0.001, 1.0);

    sollya::run_and_print_remez_f32(f, range, &poly_i, 0, "K");
}
