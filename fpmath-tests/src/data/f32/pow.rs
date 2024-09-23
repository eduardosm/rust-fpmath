use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_pow",
        gen_args,
        |(x, y)| {
            let bigx = rug::Float::with_val(24, x);
            let bigy = rug::Float::with_val(24, y);
            let tmp = rug::Float::with_val(24 * 2, rug::ops::Pow::pow(&bigx, &bigy));

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

    // x = sx * mx * 2^ex
    // log2(|x|) = log2(mx) + ex
    // ex <= log2(|x|) <= ex + 1

    // MIN <= |x|^y <= MAX
    // log2(MIN) / log2(|x|) <= y <= log2(MAX) / log2(|x|)

    let mut args = Vec::new();

    for ex in -126..=127 {
        let (min_y, max_y) = if ex == 0 {
            (-126, 127)
        } else {
            let a = 127 / ex;
            let b = -a;
            ((a.min(b) - 3).max(-126), (a.max(b) + 3).min(127))
        };

        for yi in min_y..=max_y {
            for _ in 0..100 {
                let mx = rng.gen::<u32>();
                let sx = false;
                let x = mkfloat(mx, ex, sx);

                let y = (rng.gen::<f32>() - 0.5) + (yi as f32);
                args.push((x, y));
            }
        }
    }

    for ex in -22..=-1 {
        for ey in 1..=127 {
            for _ in 0..2000 {
                let mx = rng.gen::<u32>();
                let sx = rng.gen::<bool>();
                let my = rng.gen::<u32>();
                let sy = rng.gen::<bool>();
                args.push((1.0 + mkfloat(mx, ex, sx), mkfloat(my, ey, sy)));
            }
        }
    }

    args
}
