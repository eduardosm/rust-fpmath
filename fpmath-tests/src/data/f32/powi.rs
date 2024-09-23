use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

#[derive(bincode::Encode, bincode::Decode)]
pub(crate) struct Data {
    pub(crate) x: f32,
    pub(crate) y: i32,
    pub(crate) expected: RefResult,
}

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_powi",
        gen_args,
        |(x, y)| {
            let bigx = rug::Float::with_val(24, x);
            let bigy = rug::Float::with_val(32, y);
            let tmp = rug::Float::with_val(24 * 2, rug::ops::Pow::pow(&bigx, &bigy));

            Data {
                x,
                y,
                expected: RefResult::from_rug(tmp),
            }
        },
        pb,
    );
}

fn gen_args() -> Vec<(f32, i32)> {
    let mut rng = create_prng();

    let mut args = Vec::new();

    for ex in -126..=127 {
        let (min_y, max_y) = if ex == 0 {
            (-126, 127)
        } else {
            let a = 127 / ex;
            let b = -a;
            ((a.min(b) - 3).max(-126), (a.max(b) + 3).min(127))
        };

        for y in min_y..=max_y {
            for _ in 0..100 {
                let mx = rng.gen::<u32>();
                let sx = false;
                let x = mkfloat(mx, ex, sx);

                args.push((x, i32::from(y)));
            }
        }
    }

    for ex in -22..=-1 {
        for i in (1..=31).rev() {
            for _ in 0..2000 {
                let mx = rng.gen::<u32>();
                let sx = rng.gen::<bool>();
                let x = 1.0 + mkfloat(mx, ex, sx);
                let y = ((rng.gen::<u32>() | 0x8000_0000) >> i) as i32;

                args.push((x, y));
                args.push((x, -y));
            }
        }
    }

    args
}
