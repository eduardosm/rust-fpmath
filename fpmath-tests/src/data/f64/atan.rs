use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_atan",
        gen_args,
        |x| {
            let mut tmp = dev_mpfr::Mpfr::new();
            tmp.set_prec(53 * 2);
            tmp.set_f64(x, dev_mpfr::Rnd::N);
            tmp.atan(None, dev_mpfr::Rnd::N);

            super::OneArgData {
                x,
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
        "f64_atand",
        gen_args,
        |x| {
            let mut tmp = dev_mpfr::Mpfr::new();
            tmp.set_prec(53 * 3);
            tmp.set_f64(x, dev_mpfr::Rnd::N);
            tmp.atan(None, dev_mpfr::Rnd::N);
            tmp.mul(None, Some(&conv), dev_mpfr::Rnd::N);

            super::OneArgData {
                x,
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
        "f64_atanpi",
        gen_args,
        |x| {
            let mut tmp = dev_mpfr::Mpfr::new();
            tmp.set_prec(53 * 3);
            tmp.set_f64(x, dev_mpfr::Rnd::N);
            tmp.atan(None, dev_mpfr::Rnd::N);
            tmp.div(None, Some(&conv), dev_mpfr::Rnd::N);

            super::OneArgData {
                x,
                expected: RefResult::from_mpfr(&mut tmp),
            }
        },
        pb,
    );
}

fn gen_args() -> Vec<f64> {
    let mut rng = create_prng();

    let mut args = Vec::new();
    for e in -1022..=1023 {
        args.push(mkfloat(0, e, false));
        args.push(mkfloat(0, e, true));
        args.push(mkfloat(u64::MAX, e, false));
        args.push(mkfloat(u64::MAX, e, true));

        for _ in 0..5000 {
            let m = rng.gen::<u64>();
            let s = rng.gen::<bool>();
            args.push(mkfloat(m, e, s));
        }
    }

    args
}
