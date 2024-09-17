use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

fn gen_args() -> Vec<(f64, f64)> {
    let mut rng = create_prng();

    let mut args = Vec::new();
    for ey in -1022..=1023 {
        if matches!(ey, -900..=900) && (ey & 3) != 3 {
            continue; // speed up tests
        }
        for ex in -1022..=1023 {
            if matches!(ex, -900..=900) && (ex & 3) != 3 {
                continue; // speed up tests
            }
            let my = rng.gen::<u64>();
            let mx = rng.gen::<u64>();
            args.push((mkfloat(my, ey, false), mkfloat(mx, ex, false)));
            args.push((mkfloat(my, ey, false), mkfloat(mx, ex, true)));
            args.push((mkfloat(my, ey, true), mkfloat(mx, ex, false)));
            args.push((mkfloat(my, ey, true), mkfloat(mx, ex, true)));
        }
    }

    for e in -1022..=1023 {
        for _ in 0..1000 {
            let my = rng.gen::<u64>();
            let sy = rng.gen::<bool>();
            let mx = rng.gen::<u64>();
            let sx = rng.gen::<bool>();
            args.push((mkfloat(my, e, sy), mkfloat(mx, e, sx)));

            let my = rng.gen::<u64>();
            let sy = rng.gen::<bool>();
            let mx = rng.gen::<u64>();
            let sx = rng.gen::<bool>();
            args.push((mkfloat(my, 0, sy), mkfloat(mx, e, sx)));

            let my = rng.gen::<u64>();
            let sy = rng.gen::<bool>();
            let mx = rng.gen::<u64>();
            let sx = rng.gen::<bool>();
            args.push((mkfloat(my, e, sy), mkfloat(mx, 0, sx)));
        }
    }

    args
}

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_atan2",
        gen_args,
        |(x, y)| {
            let mut bigx = dev_mpfr::Mpfr::new();
            bigx.set_prec(53);
            bigx.set_f64(x, dev_mpfr::Rnd::N);

            let mut bigy = dev_mpfr::Mpfr::new();
            bigy.set_prec(53);
            bigy.set_f64(y, dev_mpfr::Rnd::N);

            let mut tmp = dev_mpfr::Mpfr::new();
            tmp.set_prec(53 * 2);
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
    conv.set_prec(53 * 3);
    conv.const_pi(dev_mpfr::Rnd::N);
    conv.f64_div(180.0, None, dev_mpfr::Rnd::N);

    generate_data(
        "f64_atan2d",
        gen_args,
        |(x, y)| {
            let mut bigx = dev_mpfr::Mpfr::new();
            bigx.set_prec(53);
            bigx.set_f64(x, dev_mpfr::Rnd::N);

            let mut bigy = dev_mpfr::Mpfr::new();
            bigy.set_prec(53);
            bigy.set_f64(y, dev_mpfr::Rnd::N);

            let mut tmp = dev_mpfr::Mpfr::new();
            tmp.set_prec(53 * 3);
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
    conv.set_prec(53 * 3);
    conv.const_pi(dev_mpfr::Rnd::N);

    generate_data(
        "f64_atan2pi",
        gen_args,
        |(x, y)| {
            let mut bigx = dev_mpfr::Mpfr::new();
            bigx.set_prec(53);
            bigx.set_f64(x, dev_mpfr::Rnd::N);

            let mut bigy = dev_mpfr::Mpfr::new();
            bigy.set_prec(53);
            bigy.set_f64(y, dev_mpfr::Rnd::N);

            let mut tmp = dev_mpfr::Mpfr::new();
            tmp.set_prec(53 * 3);
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
