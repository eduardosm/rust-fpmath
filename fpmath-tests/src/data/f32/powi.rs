use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

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

#[derive(bincode::Encode, bincode::Decode)]
pub(crate) struct Data {
    pub(crate) x: f32,
    pub(crate) y: i32,
    pub(crate) expected: RefResult,
}

#[allow(clippy::unnecessary_fallible_conversions)]
pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_powi",
        gen_args,
        |(x, y)| {
            let mut bigx = dev_mpfr::Mpfr::new();
            bigx.set_prec(24);
            bigx.set_f32(x, dev_mpfr::Rnd::N);

            let mut bigy = dev_mpfr::Mpfr::new();
            bigy.set_prec(32);
            bigy.set_si(y.try_into().unwrap(), dev_mpfr::Rnd::N);

            let mut tmp = dev_mpfr::Mpfr::new();
            tmp.set_prec(24 * 2);
            tmp.pow(&bigx, &bigy, dev_mpfr::Rnd::N);

            Data {
                x,
                y,
                expected: RefResult::from_mpfr(&mut tmp),
            }
        },
        pb,
    );
}
