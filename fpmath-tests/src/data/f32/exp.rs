use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

fn gen_args() -> Vec<f32> {
    let mut rng = create_prng();

    let mut args = Vec::new();
    for e in -126..=9 {
        args.push(mkfloat(0, e, false));
        args.push(mkfloat(0, e, true));
        args.push(mkfloat(u32::MAX, e, false));
        args.push(mkfloat(u32::MAX, e, true));

        for _ in 0..10000 {
            let m = rng.gen::<u32>();
            let s = rng.gen::<bool>();
            args.push(mkfloat(m, e, s));
        }
    }

    for arg in -160..=160 {
        args.push(arg as f32);
    }

    args
}

#[derive(bincode::Encode, bincode::Decode)]
pub(crate) struct ExpExpM1Data {
    pub(crate) x: f32,
    pub(crate) expected_exp: RefResult,
    pub(crate) expected_expm1: RefResult,
}

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_exp",
        gen_args,
        |x| {
            let mut bigx = dev_mpfr::Mpfr::new();
            bigx.set_prec(24);
            bigx.set_f32(x, dev_mpfr::Rnd::N);

            let mut tmp_exp = dev_mpfr::Mpfr::new();
            tmp_exp.set_prec(24 * 2);
            tmp_exp.exp(Some(&bigx), dev_mpfr::Rnd::N);

            let mut tmp_expm1 = dev_mpfr::Mpfr::new();
            tmp_expm1.set_prec(24 * 2);
            tmp_expm1.expm1(Some(&bigx), dev_mpfr::Rnd::N);

            ExpExpM1Data {
                x,
                expected_exp: RefResult::from_mpfr(&mut tmp_exp),
                expected_expm1: RefResult::from_mpfr(&mut tmp_expm1),
            }
        },
        pb,
    );
}

pub(crate) fn gen_data_2(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_exp2",
        gen_args,
        |x| {
            let mut tmp = dev_mpfr::Mpfr::new();
            tmp.set_prec(24 * 2);
            tmp.set_f32(x, dev_mpfr::Rnd::N);
            tmp.exp2(None, dev_mpfr::Rnd::N);

            super::OneArgData {
                x,
                expected: RefResult::from_mpfr(&mut tmp),
            }
        },
        pb,
    );
}

pub(crate) fn gen_data_10(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_exp10",
        gen_args,
        |x| {
            let mut tmp = dev_mpfr::Mpfr::new();
            tmp.set_prec(24 * 2);
            tmp.set_f32(x, dev_mpfr::Rnd::N);
            tmp.exp10(None, dev_mpfr::Rnd::N);

            super::OneArgData {
                x,
                expected: RefResult::from_mpfr(&mut tmp),
            }
        },
        pb,
    );
}
