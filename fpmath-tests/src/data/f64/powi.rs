use rand::Rng as _;

use super::{mkfloat, RefResult, RUG_PREC};
use crate::data::{create_prng, generate_data};

#[derive(bincode::Encode, bincode::Decode)]
pub(crate) struct Data {
    pub(crate) x: f64,
    pub(crate) y: i32,
    pub(crate) expected: RefResult,
}

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_powi",
        gen_args,
        |(x, y)| {
            let bigx = rug::Float::with_val(53, x);
            let bigy = rug::Float::with_val(32, y);
            let tmp = rug::Float::with_val(RUG_PREC, rug::ops::Pow::pow(&bigx, &bigy));

            Data {
                x,
                y,
                expected: RefResult::from_rug(tmp),
            }
        },
        pb,
    );
}

fn gen_args() -> Vec<(f64, i32)> {
    let mut rng = create_prng();

    let mut args = Vec::new();

    for ex in -1022..=1023 {
        let (min_y, max_y) = if ex == 0 {
            (-1022, 1023)
        } else {
            let a = 1023 / ex;
            let b = -a;
            ((a.min(b) - 3).max(-1022), (a.max(b) + 3).min(1023))
        };

        for y in min_y..=max_y {
            for _ in 0..100 {
                let mx = rng.random::<u64>();
                let sx = false;
                let x = mkfloat(mx, ex, sx);

                args.push((x, i32::from(y)));
            }
        }
    }

    for ex in -51..=-1 {
        for i in (1..=31).rev() {
            for _ in 0..1000 {
                let mx = rng.random::<u64>();
                let sx = rng.random::<bool>();
                let x = 1.0 + mkfloat(mx, ex, sx);
                let y = ((rng.random::<u32>() | 0x8000_0000) >> i) as i32;

                args.push((x, y));
                args.push((x, -y));
            }
        }
    }

    args
}
