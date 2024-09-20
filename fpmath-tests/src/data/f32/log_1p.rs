use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_log_1p",
        gen_args,
        |x| {
            let mut tmp = dev_mpfr::Mpfr::new();
            tmp.set_prec(24 * 2);
            tmp.set_f32(x, dev_mpfr::Rnd::N);
            tmp.log1p(None, dev_mpfr::Rnd::N);

            super::OneArgData {
                x,
                expected: RefResult::from_mpfr(&mut tmp),
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

    // < 0
    for e in -126..=-1 {
        args.push(mkfloat(0, e, true));
        args.push(mkfloat(u32::MAX, e, true));

        for _ in 0..1000 {
            let m = rng.gen::<u32>();
            args.push(mkfloat(m, e, true));
        }
    }

    // subnormals
    for i in 0..23 {
        args.push(f32::from_bits(1 << i));
        args.push(f32::from_bits((1 << (i + 1)) - 1));
    }

    args
}
