use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

fn gen_args() -> Vec<f64> {
    let mut rng = create_prng();

    let mut args = Vec::new();
    for e in -1022..=-1 {
        args.push(mkfloat(0, e, false));
        args.push(mkfloat(0, e, true));
        args.push(mkfloat(u64::MAX, e, false));
        args.push(mkfloat(u64::MAX, e, true));

        for _ in 0..10000 {
            let m = rng.gen::<u64>();
            let s = rng.gen::<bool>();
            args.push(mkfloat(m, e, s));
        }
    }

    args
}

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_atanh",
        gen_args,
        |x| {
            let mut tmp = dev_mpfr::Mpfr::new();
            tmp.set_prec(53 * 2);
            tmp.set_f64(x, dev_mpfr::Rnd::N);
            tmp.atanh(None, dev_mpfr::Rnd::N);

            super::OneArgData {
                x,
                expected: RefResult::from_mpfr(&mut tmp),
            }
        },
        pb,
    );
}
