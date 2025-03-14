use rand::Rng as _;

use super::{mkfloat, RefResult, RUG_PREC};
use crate::data::{create_prng, generate_data};

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_atan2",
        gen_args,
        |(y, x)| {
            let bigy = rug::Float::with_val(53, y);
            let bigx = rug::Float::with_val(53, x);
            let tmp = rug::Float::with_val(RUG_PREC, bigy.atan2_ref(&bigx));

            super::TwoArgData {
                x,
                y,
                expected: RefResult::from_rug(tmp),
            }
        },
        pb,
    );
}

pub(crate) fn gen_data_d(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_atan2d",
        gen_args,
        |(y, x)| {
            let bigy = rug::Float::with_val(53, y);
            let bigx = rug::Float::with_val(53, x);
            let tmp = rug::Float::with_val(RUG_PREC, bigy.atan2_u_ref(&bigx, 360));

            super::TwoArgData {
                x,
                y,
                expected: RefResult::from_rug(tmp),
            }
        },
        pb,
    );
}

pub(crate) fn gen_data_pi(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_atan2pi",
        gen_args,
        |(y, x)| {
            let bigy = rug::Float::with_val(53, y);
            let bigx = rug::Float::with_val(53, x);
            let tmp = rug::Float::with_val(RUG_PREC, bigy.atan2_pi_ref(&bigx));

            super::TwoArgData {
                x,
                y,
                expected: RefResult::from_rug(tmp),
            }
        },
        pb,
    );
}

fn gen_args() -> Vec<(f64, f64)> {
    let mut rng = create_prng();

    let mut args = Vec::new();
    for ey in -1022..=1023 {
        if matches!(ey, -900..=900) && (ey & 3) != 3 {
            continue; // speed up tests
        }
        for ex in -1022..=1023 {
            if matches!(ex, -900..=900) && (ex & 3) != 3 {
                continue; // speed up tests
            }
            let my = rng.random::<u64>();
            let mx = rng.random::<u64>();
            args.push((mkfloat(my, ey, false), mkfloat(mx, ex, false)));
            args.push((mkfloat(my, ey, false), mkfloat(mx, ex, true)));
        }
    }

    for e in -1022..=1023 {
        for _ in 0..1000 {
            let my = rng.random::<u64>();
            let mx = rng.random::<u64>();
            let sx = rng.random::<bool>();
            args.push((mkfloat(my, e, false), mkfloat(mx, e, sx)));

            let my = rng.random::<u64>();
            let mx = rng.random::<u64>();
            let sx = rng.random::<bool>();
            args.push((mkfloat(my, 0, false), mkfloat(mx, e, sx)));

            let my = rng.random::<u64>();
            let mx = rng.random::<u64>();
            let sx = rng.random::<bool>();
            args.push((mkfloat(my, e, false), mkfloat(mx, 0, sx)));
        }
    }

    args
}
