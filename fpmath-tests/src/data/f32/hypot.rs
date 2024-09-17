use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

fn gen_args() -> Vec<(f32, f32)> {
    let mut rng = create_prng();

    let mut args = Vec::new();

    for ex in -126..=127 {
        for ey in -126..=127 {
            for _ in 0..5 {
                let mx = rng.gen::<u32>();
                let sx = rng.gen::<bool>();
                let my = rng.gen::<u32>();
                let sy = rng.gen::<bool>();
                args.push((mkfloat(mx, ex, sx), mkfloat(my, ey, sy)));
            }
        }
    }

    for e in -126..=127 {
        for _ in 0..5000 {
            let mx = rng.gen::<u32>();
            let sx = rng.gen::<bool>();
            let my = rng.gen::<u32>();
            let sy = rng.gen::<bool>();
            args.push((mkfloat(mx, e, sx), mkfloat(my, e, sy)));

            let mx = rng.gen::<u32>();
            let sx = rng.gen::<bool>();
            let my = rng.gen::<u32>();
            let sy = rng.gen::<bool>();
            args.push((mkfloat(mx, 0, sx), mkfloat(my, e, sy)));

            let mx = rng.gen::<u32>();
            let sx = rng.gen::<bool>();
            let my = rng.gen::<u32>();
            let sy = rng.gen::<bool>();
            args.push((mkfloat(mx, e, sx), mkfloat(my, 0, sy)));
        }
    }

    args
}

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_hypot",
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
            tmp.hypot(&bigx, &bigy, dev_mpfr::Rnd::N);

            super::TwoArgData {
                x,
                y,
                expected: RefResult::from_mpfr(&mut tmp),
            }
        },
        pb,
    );
}
