use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_log_1p",
        gen_args,
        |x| {
            let tmp = rug::Float::with_val(53 * 2, x).ln_1p();

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
    for e in -100..=100 {
        args.push(mkfloat(0, e, false));
        args.push(mkfloat(u64::MAX, e, false));

        for _ in 0..9000 {
            let m = rng.random::<u64>();
            args.push(mkfloat(m, e, false));
        }
    }
    for e in -1022..=1023 {
        for _ in 0..1000 {
            let m = rng.random::<u64>();
            args.push(mkfloat(m, e, false));
        }
    }

    for arg in 1..=10000 {
        args.push(arg as f64);
    }

    args.push(f64::MIN_POSITIVE);
    args.push(f64::MAX);

    // 1 < x < 0
    for e in -1022..=-1 {
        args.push(mkfloat(0, e, true));
        args.push(mkfloat(u64::MAX, e, true));

        for _ in 0..1000 {
            let m = rng.random::<u64>();
            args.push(mkfloat(m, e, true));
        }
    }

    // subnormals
    for i in 0..52 {
        args.push(f64::from_bits(1 << i));
        args.push(-f64::from_bits(1 << i));
        args.push(f64::from_bits((1 << (i + 1)) - 1));
        args.push(-f64::from_bits((1 << (i + 1)) - 1));
    }

    args
}
