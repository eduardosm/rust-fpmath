use super::{sinh_cosh::gen_args, RefResult};
use crate::data::generate_data;

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_tanh",
        gen_args,
        |x| {
            let tmp = rug::Float::with_val(53 * 2, x).tanh();

            super::OneArgData {
                x,
                expected: RefResult::from_rug(tmp),
            }
        },
        pb,
    );
}
