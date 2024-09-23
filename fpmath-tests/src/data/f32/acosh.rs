use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_acosh",
        gen_args,
        |x| {
            let tmp = rug::Float::with_val(24 * 2, x).acosh();

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
    for e in 0..=127 {
        args.push(mkfloat(0, e, false));
        args.push(mkfloat(u32::MAX, e, false));

        for _ in 0..10000 {
            let m = rng.gen::<u32>();
            args.push(mkfloat(m, e, false));
        }
    }

    args
}
