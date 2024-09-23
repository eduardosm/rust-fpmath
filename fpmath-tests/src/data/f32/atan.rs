use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_atan",
        gen_args,
        |x| {
            let tmp = rug::Float::with_val(24 * 2, x).atan();

            super::OneArgData {
                x,
                expected: RefResult::from_rug(tmp),
            }
        },
        pb,
    );
}

pub(crate) fn gen_data_d(pb: indicatif::ProgressBar) {
    let conv = 180u8 / rug::Float::with_val(24 * 3, rug::float::Constant::Pi);

    generate_data(
        "f32_atand",
        gen_args,
        |x| {
            let tmp = rug::Float::with_val(24 * 3, x).atan() * &conv;

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
        "f32_atanpi",
        gen_args,
        |x| {
            let tmp = rug::Float::with_val(24 * 2, x).atan_pi();

            super::OneArgData {
                x,
                expected: RefResult::from_rug(tmp),
            }
        },
        pb,
    );
}

fn gen_args() -> Vec<f32> {
    let mut rng = create_prng();

    let mut args = Vec::new();
    for e in -126..=127 {
        args.push(mkfloat(0, e, false));
        args.push(mkfloat(0, e, true));
        args.push(mkfloat(u32::MAX, e, false));
        args.push(mkfloat(u32::MAX, e, true));

        for _ in 0..5000 {
            let m = rng.gen::<u32>();
            let s = rng.gen::<bool>();
            args.push(mkfloat(m, e, s));
        }
    }

    args
}
