use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

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

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_pow",
        gen_args,
        |(x, y)| {
            let mut bigx = dev_mpfr::Mpfr::new();
            bigx.set_prec(24);
            bigx.set_f32(x, dev_mpfr::Rnd::N);

            let mut bigy = dev_mpfr::Mpfr::new();
            bigy.set_prec(24);
            bigy.set_f32(y, dev_mpfr::Rnd::N);

            let mut tmp = dev_mpfr::Mpfr::new();
            tmp.set_prec(24 * 2);
            tmp.pow(&bigx, &bigy, dev_mpfr::Rnd::N);

            super::TwoArgData {
                x,
                y,
                expected: RefResult::from_mpfr(&mut tmp),
            }
        },
        pb,
    );
}
