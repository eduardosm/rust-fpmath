use rand::Rng as _;

use super::{mkfloat, RefResult, RUG_PREC};
use crate::data::{create_prng, generate_data};

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_sind_cosd",
        gen_args,
        |x| {
            let bigx = rug::Float::with_val(53, x);
            let tmp_sin = rug::Float::with_val(RUG_PREC, bigx.sin_u_ref(360));
            let tmp_cos = rug::Float::with_val(RUG_PREC, bigx.cos_u_ref(360));

            super::SinCosData {
                x,
                expected_sin: RefResult::from_rug(tmp_sin),
                expected_cos: RefResult::from_rug(tmp_cos),
            }
        },
        pb,
    );
}

pub(super) fn gen_args() -> Vec<f64> {
    let mut rng = create_prng();

    let mut args = Vec::new();
    for e in -1022..=1023 {
        args.push(mkfloat(0, e, false));
        args.push(mkfloat(u64::MAX, e, false));

        for _ in 0..1000 {
            let m = rng.random::<u64>();
            args.push(mkfloat(m, e, false));
        }
    }

    for arg in 1..=20000 {
        args.push((arg as f64) * 0.5);
    }

    args
}
