use super::{sin_cos::gen_args, RefResult};
use crate::data::generate_data;

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_tan",
        gen_args,
        |x| {
            let tmp = rug::Float::with_val(24 * 2, x).tan();

            super::OneArgData {
                x,
                expected: RefResult::from_rug(tmp),
            }
        },
        pb,
    );
}
