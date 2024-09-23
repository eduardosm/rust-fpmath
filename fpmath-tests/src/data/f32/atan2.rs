use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_atan2",
        gen_args,
        |(x, y)| {
            let bigx = rug::Float::with_val(24, x);
            let bigy = rug::Float::with_val(24, y);
            let tmp = rug::Float::with_val(24 * 2, bigx.atan2_ref(&bigy));

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
    let conv = 180u8 / rug::Float::with_val(24 * 3, rug::float::Constant::Pi);

    generate_data(
        "f32_atan2d",
        gen_args,
        |(x, y)| {
            let bigx = rug::Float::with_val(24, x);
            let bigy = rug::Float::with_val(24, y);
            let tmp = rug::Float::with_val(24 * 3, bigx.atan2_ref(&bigy)) * &conv;

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
        "f32_atan2pi",
        gen_args,
        |(x, y)| {
            let bigx = rug::Float::with_val(24, x);
            let bigy = rug::Float::with_val(24, y);
            let tmp = rug::Float::with_val(24 * 2, bigx.atan2_pi_ref(&bigy));

            super::TwoArgData {
                x,
                y,
                expected: RefResult::from_rug(tmp),
            }
        },
        pb,
    );
}

fn gen_args() -> Vec<(f32, f32)> {
    let mut rng = create_prng();

    let mut args = Vec::new();
    for ey in -126..=127 {
        for ex in -126..=127 {
            let my = rng.gen::<u32>();
            let mx = rng.gen::<u32>();
            args.push((mkfloat(my, ey, false), mkfloat(mx, ex, false)));
            args.push((mkfloat(my, ey, false), mkfloat(mx, ex, true)));
            args.push((mkfloat(my, ey, true), mkfloat(mx, ex, false)));
            args.push((mkfloat(my, ey, true), mkfloat(mx, ex, true)));
        }
    }

    for e in -126..=127 {
        for _ in 0..5000 {
            let my = rng.gen::<u32>();
            let sy = rng.gen::<bool>();
            let mx = rng.gen::<u32>();
            let sx = rng.gen::<bool>();
            args.push((mkfloat(my, e, sy), mkfloat(mx, e, sx)));

            let my = rng.gen::<u32>();
            let sy = rng.gen::<bool>();
            let mx = rng.gen::<u32>();
            let sx = rng.gen::<bool>();
            args.push((mkfloat(my, 0, sy), mkfloat(mx, e, sx)));

            let my = rng.gen::<u32>();
            let sy = rng.gen::<bool>();
            let mx = rng.gen::<u32>();
            let sx = rng.gen::<bool>();
            args.push((mkfloat(my, e, sy), mkfloat(mx, 0, sx)));
        }
    }

    args
}
