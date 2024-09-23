use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_atan",
        gen_args,
        |x| {
            let tmp = rug::Float::with_val(53 * 2, x).atan();

            super::OneArgData {
                x,
                expected: RefResult::from_rug(tmp),
            }
        },
        pb,
    );
}

pub(crate) fn gen_data_d(pb: indicatif::ProgressBar) {
    let conv = 180u8 / rug::Float::with_val(53 * 3, rug::float::Constant::Pi);

    generate_data(
        "f64_atand",
        gen_args,
        |x| {
            let tmp = rug::Float::with_val(53 * 3, x).atan() * &conv;

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
            let tmp = rug::Float::with_val(53 * 2, x).atan_pi();

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
        args.push(mkfloat(0, e, true));
        args.push(mkfloat(u64::MAX, e, false));
        args.push(mkfloat(u64::MAX, e, true));

        for _ in 0..5000 {
            let m = rng.gen::<u64>();
            let s = rng.gen::<bool>();
            args.push(mkfloat(m, e, s));
        }
    }

    args
}
