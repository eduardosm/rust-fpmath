use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_pow",
        gen_args,
        |(x, y)| {
            let bigx = rug::Float::with_val(53, x);
            let bigy = rug::Float::with_val(53, y);
            let tmp = rug::Float::with_val(53 * 2, rug::ops::Pow::pow(&bigx, &bigy));

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

    // x = sx * mx * 2^ex
    // log2(|x|) = log2(mx) + ex
    // ex <= log2(|x|) <= ex + 1

    // MIN <= |x|^y <= MAX
    // log2(MIN) / log2(|x|) <= y <= log2(MAX) / log2(|x|)

    let mut args = Vec::new();

    for ex in -1022..=1023 {
        let (min_y, max_y) = if ex == 0 {
            (-1022, 1023)
        } else {
            let a = 1023 / ex;
            let b = -a;
            ((a.min(b) - 3).max(-1022), (a.max(b) + 3).min(1023))
        };

        for yi in min_y..=max_y {
            for _ in 0..50 {
                let mx = rng.random::<u64>();
                let sx = false;
                let x = mkfloat(mx, ex, sx);

                let y = (rng.random::<f64>() - 0.5) + (yi as f64);
                args.push((x, y));
            }
        }
    }

    for ex in -51..=-1 {
        for ey in 1..=1023 {
            for _ in 0..100 {
                let mx = rng.random::<u64>();
                let sx = rng.random::<bool>();
                let my = rng.random::<u64>();
                let sy = rng.random::<bool>();
                args.push((1.0 + mkfloat(mx, ex, sx), mkfloat(my, ey, sy)));
            }
        }
    }

    args
}
