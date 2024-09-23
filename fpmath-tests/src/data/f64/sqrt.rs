use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_sqrt",
        gen_args,
        |x| {
            let tmp = rug::Float::with_val(53 * 2, x).sqrt();

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

    for e in -200..=200 {
        args.push(mkfloat(0, e, false));
        args.push(mkfloat(u64::MAX, e, false));

        for _ in 0..5000 {
            let m = rng.gen::<u64>();
            args.push(mkfloat(m, e, false));
        }
    }
    for e in -1022..=1023 {
        for _ in 0..5000 {
            let m = rng.gen::<u64>();
            args.push(mkfloat(m, e, false));
        }
    }

    for arg in 1..=10000 {
        args.push(arg as f64);
    }

    args.push(f64::MIN_POSITIVE);
    args.push(f64::MAX);

    // subnormals
    for i in 0..52 {
        args.push(f64::from_bits(1 << i));
        args.push(f64::from_bits((1 << (i + 1)) - 1));
    }

    args
}
