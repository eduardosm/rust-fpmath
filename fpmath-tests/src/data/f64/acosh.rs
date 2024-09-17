use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

fn gen_args() -> Vec<f64> {
    let mut rng = create_prng();

    let mut args = Vec::new();
    for e in 0..=100 {
        args.push(mkfloat(0, e, false));
        args.push(mkfloat(u64::MAX, e, false));

        for _ in 0..9000 {
            let m = rng.gen::<u64>();
            args.push(mkfloat(m, e, false));
        }
    }

    for e in 0..=1023 {
        args.push(mkfloat(0, e, false));
        args.push(mkfloat(u64::MAX, e, false));

        for _ in 0..1000 {
            let m = rng.gen::<u64>();
            args.push(mkfloat(m, e, false));
        }
    }

    args
}

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_acosh",
        gen_args,
        |x| {
            let mut tmp = dev_mpfr::Mpfr::new();
            tmp.set_prec(53 * 2);
            tmp.set_f64(x, dev_mpfr::Rnd::N);
            tmp.acosh(None, dev_mpfr::Rnd::N);

            super::OneArgData {
                x,
                expected: RefResult::from_mpfr(&mut tmp),
            }
        },
        pb,
    );
}
