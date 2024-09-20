use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_atan2",
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
            tmp.atan2(&bigx, &bigy, dev_mpfr::Rnd::N);

            super::TwoArgData {
                x,
                y,
                expected: RefResult::from_mpfr(&mut tmp),
            }
        },
        pb,
    );
}

pub(crate) fn gen_data_d(pb: indicatif::ProgressBar) {
    let mut conv = dev_mpfr::Mpfr::new();
    conv.set_prec(24 * 3);
    conv.const_pi(dev_mpfr::Rnd::N);
    conv.f64_div(180.0, None, dev_mpfr::Rnd::N);

    generate_data(
        "f32_atan2d",
        gen_args,
        |(x, y)| {
            let mut bigx = dev_mpfr::Mpfr::new();
            bigx.set_prec(24);
            bigx.set_f32(x, dev_mpfr::Rnd::N);

            let mut bigy = dev_mpfr::Mpfr::new();
            bigy.set_prec(24);
            bigy.set_f32(y, dev_mpfr::Rnd::N);

            let mut tmp = dev_mpfr::Mpfr::new();
            tmp.set_prec(24 * 3);
            tmp.atan2(&bigx, &bigy, dev_mpfr::Rnd::N);
            tmp.mul(None, Some(&conv), dev_mpfr::Rnd::N);

            super::TwoArgData {
                x,
                y,
                expected: RefResult::from_mpfr(&mut tmp),
            }
        },
        pb,
    );
}

pub(crate) fn gen_data_pi(pb: indicatif::ProgressBar) {
    let mut conv = dev_mpfr::Mpfr::new();
    conv.set_prec(24 * 3);
    conv.const_pi(dev_mpfr::Rnd::N);

    generate_data(
        "f32_atan2pi",
        gen_args,
        |(x, y)| {
            let mut bigx = dev_mpfr::Mpfr::new();
            bigx.set_prec(24);
            bigx.set_f32(x, dev_mpfr::Rnd::N);

            let mut bigy = dev_mpfr::Mpfr::new();
            bigy.set_prec(24);
            bigy.set_f32(y, dev_mpfr::Rnd::N);

            let mut tmp = dev_mpfr::Mpfr::new();
            tmp.set_prec(24 * 3);
            tmp.atan2(&bigx, &bigy, dev_mpfr::Rnd::N);
            tmp.div(None, Some(&conv), dev_mpfr::Rnd::N);

            super::TwoArgData {
                x,
                y,
                expected: RefResult::from_mpfr(&mut tmp),
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
