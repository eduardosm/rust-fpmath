use super::{sinh_cosh::gen_args, RefResult, RUG_PREC};
use crate::data::generate_data;

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_tanh",
        gen_args,
        |x| {
            let tmp = rug::Float::with_val(RUG_PREC, x).tanh();

            super::OneArgData {
                x,
                expected: RefResult::from_rug(tmp),
            }
        },
        pb,
    );
}
