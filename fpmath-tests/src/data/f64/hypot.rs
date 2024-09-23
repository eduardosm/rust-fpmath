use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_hypot",
        gen_args,
        |(x, y)| {
            let bigx = rug::Float::with_val(53, x);
            let bigy = rug::Float::with_val(53, y);
            let tmp = rug::Float::with_val(53 * 2, bigx.hypot_ref(&bigy));

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

    for ex in -1022..=1023 {
        if matches!(ex, -900..=900) && (ex & 3) != 3 {
            continue; // speed up tests
        }
        for ey in -1022..=1023 {
            for _ in 0..5 {
                let mx = rng.gen::<u64>();
                let sx = rng.gen::<bool>();
                let my = rng.gen::<u64>();
                let sy = rng.gen::<bool>();
                args.push((mkfloat(mx, ex, sx), mkfloat(my, ey, sy)));
            }
        }
    }

    for e in -1022..=1023 {
        for _ in 0..1000 {
            let mx = rng.gen::<u64>();
            let sx = rng.gen::<bool>();
            let my = rng.gen::<u64>();
            let sy = rng.gen::<bool>();
            args.push((mkfloat(mx, e, sx), mkfloat(my, e, sy)));

            let mx = rng.gen::<u64>();
            let sx = rng.gen::<bool>();
            let my = rng.gen::<u64>();
            let sy = rng.gen::<bool>();
            args.push((mkfloat(mx, 0, sx), mkfloat(my, e, sy)));

            let mx = rng.gen::<u64>();
            let sx = rng.gen::<bool>();
            let my = rng.gen::<u64>();
            let sy = rng.gen::<bool>();
            args.push((mkfloat(mx, e, sx), mkfloat(my, 0, sy)));
        }
    }

    args
}
