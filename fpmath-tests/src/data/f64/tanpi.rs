use super::{sind_cosd::gen_args, RefResult};
use crate::data::generate_data;

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    let mut k_2 = dev_mpfr::Mpfr::new();
    k_2.set_prec(128);
    k_2.set_ui(2, dev_mpfr::Rnd::N);

    let mut conv = dev_mpfr::Mpfr::new();
    conv.set_prec(2048);
    conv.const_pi(dev_mpfr::Rnd::N);

    generate_data(
        "f64_tanpi",
        gen_args,
        |x| {
            let mut tmp_arg = dev_mpfr::Mpfr::new();
            tmp_arg.set_prec(2048);
            tmp_arg.set_f64(x, dev_mpfr::Rnd::N);
            tmp_arg.fmod(None, Some(&k_2), dev_mpfr::Rnd::N);
            tmp_arg.mul(Some(&conv), None, dev_mpfr::Rnd::N);

            let mut tmp_tan = dev_mpfr::Mpfr::new();
            tmp_tan.set_prec(53 * 2);
            tmp_tan.set_f64(x, dev_mpfr::Rnd::N);
            tmp_tan.tan(Some(&tmp_arg), dev_mpfr::Rnd::N);

            if tmp_tan.get_f64(dev_mpfr::Rnd::N).abs() > f64::MAX {
                tmp_tan.set_inf(x.is_sign_negative());
            }

            super::OneArgData {
                x,
                expected: RefResult::from_mpfr(&mut tmp_tan),
            }
        },
        pb,
    );
}
