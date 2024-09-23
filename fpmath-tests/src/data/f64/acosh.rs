use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_acosh",
        gen_args,
        |x| {
            let tmp = rug::Float::with_val(53 * 2, x).acosh();

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
