use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_sinh_cosh",
        gen_args,
        |x| {
            let (tmp_sinh, tmp_cosh) =
                rug::Float::with_val(53 * 2, x).sinh_cosh(rug::Float::new(53 * 2));

            super::SinCosData {
                x,
                expected_sin: RefResult::from_rug(tmp_sinh),
                expected_cos: RefResult::from_rug(tmp_cosh),
            }
        },
        pb,
    );
}

pub(super) fn gen_args() -> Vec<f64> {
    let mut rng = create_prng();

    let mut args = Vec::new();
    for e in -1022..-200 {
        args.push(mkfloat(0, e, false));
        args.push(mkfloat(u64::MAX, e, false));

        for _ in 0..100 {
            let m = rng.random::<u64>();
            args.push(mkfloat(m, e, false));
        }
    }
    for e in -200..=12 {
        args.push(mkfloat(0, e, false));
        args.push(mkfloat(u64::MAX, e, false));

        for _ in 0..1000 {
            let m = rng.random::<u64>();
            args.push(mkfloat(m, e, false));
        }
    }

    for arg in 1..=1000 {
        args.push(arg as f64);
    }

    args
}
