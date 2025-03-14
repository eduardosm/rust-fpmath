use rand::Rng as _;

use super::{mkfloat, RefResult, RUG_PREC};
use crate::data::{create_prng, generate_data};

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_atan",
        gen_args,
        |x| {
            let tmp = rug::Float::with_val(RUG_PREC, x).atan();

            super::OneArgData {
                x,
                expected: RefResult::from_rug(tmp),
            }
        },
        pb,
    );
}

pub(crate) fn gen_data_d(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_atand",
        gen_args,
        |x| {
            let tmp = rug::Float::with_val(RUG_PREC, x).atan_u(360);

            super::OneArgData {
                x,
                expected: RefResult::from_rug(tmp),
            }
        },
        pb,
    );
}

pub(crate) fn gen_data_pi(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_atanpi",
        gen_args,
        |x| {
            let tmp = rug::Float::with_val(RUG_PREC, x).atan_pi();

            super::OneArgData {
                x,
                expected: RefResult::from_rug(tmp),
            }
        },
        pb,
    );
}

fn gen_args() -> Vec<f64> {
    let mut rng = create_prng();

    let mut args = Vec::new();
    for e in -1022..=1023 {
        args.push(mkfloat(0, e, false));
        args.push(mkfloat(u64::MAX, e, false));

        for _ in 0..5000 {
            let m = rng.random::<u64>();
            args.push(mkfloat(m, e, false));
        }
    }

    args
}
