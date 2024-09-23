use super::{sind_cosd::gen_args, RefResult};
use crate::data::generate_data;

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_tanpi",
        gen_args,
        |x| {
            let mut tmp = rug::Float::with_val(24 * 2, x).tan_pi();

            if !(f32::MIN..=f32::MAX).contains(&tmp) {
                rug::Assign::assign(
                    &mut tmp,
                    if x.is_sign_positive() {
                        rug::float::Special::Infinity
                    } else {
                        rug::float::Special::NegInfinity
                    },
                );
            }

            super::OneArgData {
                x,
                expected: RefResult::from_rug(tmp),
            }
        },
        pb,
    );
}
