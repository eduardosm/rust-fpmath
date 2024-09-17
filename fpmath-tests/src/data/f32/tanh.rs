use super::{sinh_cosh::gen_args, RefResult};
use crate::data::generate_data;

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_tanh",
        gen_args,
        |x| {
            let mut tmp = dev_mpfr::Mpfr::new();
            tmp.set_prec(24 * 2);
            tmp.set_f32(x, dev_mpfr::Rnd::N);
            tmp.tanh(None, dev_mpfr::Rnd::N);

            super::OneArgData {
                x,
                expected: RefResult::from_mpfr(&mut tmp),
            }
        },
        pb,
    );
}
