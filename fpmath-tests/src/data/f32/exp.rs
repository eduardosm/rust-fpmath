use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

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
            let bigx = rug::Float::with_val(24, x);
            let tmp_exp = rug::Float::with_val(24 * 2, bigx.exp_ref());
            let tmp_expm1 = rug::Float::with_val(24 * 2, bigx.exp_m1_ref());

            ExpExpM1Data {
                x,
                expected_exp: RefResult::from_rug(tmp_exp),
                expected_expm1: RefResult::from_rug(tmp_expm1),
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
            let tmp = rug::Float::with_val(24 * 2, x).exp2();

            super::OneArgData {
                x,
                expected: RefResult::from_rug(tmp),
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
            let tmp = rug::Float::with_val(24 * 2, x).exp10();

            super::OneArgData {
                x,
                expected: RefResult::from_rug(tmp),
            }
        },
        pb,
    );
}

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
