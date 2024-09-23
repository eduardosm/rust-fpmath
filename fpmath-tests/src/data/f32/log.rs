use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_log",
        gen_args,
        |x| {
            let tmp = rug::Float::with_val(24 * 2, x).ln();

            super::OneArgData {
                x,
                expected: RefResult::from_rug(tmp),
            }
        },
        pb,
    );
}

pub(crate) fn gen_data_2(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_log2",
        gen_args,
        |x| {
            let tmp = rug::Float::with_val(24 * 2, x).log2();

            super::OneArgData {
                x,
                expected: RefResult::from_rug(tmp),
            }
        },
        pb,
    );
}

pub(crate) fn gen_data_10(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_log10",
        gen_args,
        |x| {
            let tmp = rug::Float::with_val(24 * 2, x).log10();

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
        args.push(mkfloat(u32::MAX, e, false));

        for _ in 0..10000 {
            let m = rng.gen::<u32>();
            args.push(mkfloat(m, e, false));
        }
    }

    for arg in 1..=10000 {
        args.push(arg as f32);
    }

    args.push(f32::MIN_POSITIVE);
    args.push(f32::MAX);

    // subnormals
    for i in 0..23 {
        args.push(f32::from_bits(1 << i));
        args.push(f32::from_bits((1 << (i + 1)) - 1));
    }

    args
}
