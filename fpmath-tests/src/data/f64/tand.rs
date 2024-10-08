use super::{sind_cosd::gen_args, RefResult};
use crate::data::generate_data;

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_tand",
        gen_args,
        |x| {
            let tmp_arg = (rug::Float::with_val(53 * 3, x) % 360u16) / 180u8;
            let mut tmp_tan = rug::Float::with_val(53 * 2, tmp_arg.tan_pi());

            if !(f64::MIN..=f64::MAX).contains(&tmp_tan) {
                rug::Assign::assign(
                    &mut tmp_tan,
                    if x.is_sign_positive() {
                        rug::float::Special::Infinity
                    } else {
                        rug::float::Special::NegInfinity
                    },
                );
            }

            super::OneArgData {
                x,
                expected: RefResult::from_rug(tmp_tan),
            }
        },
        pb,
    );
}
