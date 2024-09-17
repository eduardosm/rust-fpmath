use super::{sind_cosd::gen_args, RefResult};
use crate::data::generate_data;

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    let mut k_360 = dev_mpfr::Mpfr::new();
    k_360.set_prec(64);
    k_360.set_ui(360, dev_mpfr::Rnd::N);

    let mut conv = dev_mpfr::Mpfr::new();
    conv.set_prec(512);
    conv.const_pi(dev_mpfr::Rnd::N);
    conv.div_f64(None, 180.0, dev_mpfr::Rnd::N);

    generate_data(
        "f32_tand",
        gen_args,
        |x| {
            let mut tmp_arg = dev_mpfr::Mpfr::new();
            tmp_arg.set_prec(512);
            tmp_arg.set_f32(x, dev_mpfr::Rnd::N);
            tmp_arg.fmod(None, Some(&k_360), dev_mpfr::Rnd::N);
            tmp_arg.mul(Some(&conv), None, dev_mpfr::Rnd::N);

            let mut tmp_tan = dev_mpfr::Mpfr::new();
            tmp_tan.set_prec(24 * 2);
            tmp_tan.set_f32(x, dev_mpfr::Rnd::N);
            tmp_tan.tan(Some(&tmp_arg), dev_mpfr::Rnd::N);

            if tmp_tan.get_f32(dev_mpfr::Rnd::N).abs() > f32::MAX {
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
