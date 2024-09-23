use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

#[derive(bincode::Encode, bincode::Decode)]
pub(crate) struct Data {
    pub(crate) x: f32,
    pub(crate) expected_asin: RefResult,
    pub(crate) expected_acos: RefResult,
}

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_asin_acos",
        gen_args,
        |x| {
            let bigx = rug::Float::with_val(24, x);
            let tmp_asin = rug::Float::with_val(24 * 2, bigx.asin_ref());
            let tmp_acos = rug::Float::with_val(24 * 2, bigx.acos_ref());

            Data {
                x,
                expected_asin: RefResult::from_rug(tmp_asin),
                expected_acos: RefResult::from_rug(tmp_acos),
            }
        },
        pb,
    );
}

pub(crate) fn gen_data_d(pb: indicatif::ProgressBar) {
    let conv = 180u8 / rug::Float::with_val(24 * 3, rug::float::Constant::Pi);

    generate_data(
        "f32_asind_acosd",
        gen_args,
        |x| {
            let bigx = rug::Float::with_val(24, x);
            let tmp_asin = rug::Float::with_val(24 * 3, bigx.asin_ref()) * &conv;
            let tmp_acos = rug::Float::with_val(24 * 3, bigx.acos_ref()) * &conv;

            Data {
                x,
                expected_asin: RefResult::from_rug(tmp_asin),
                expected_acos: RefResult::from_rug(tmp_acos),
            }
        },
        pb,
    );
}

pub(crate) fn gen_data_pi(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_asinpi_acospi",
        gen_args,
        |x| {
            let bigx = rug::Float::with_val(24, x);
            let tmp_asin = rug::Float::with_val(24 * 2, bigx.asin_pi_ref());
            let tmp_acos = rug::Float::with_val(24 * 2, bigx.acos_pi_ref());

            Data {
                x,
                expected_asin: RefResult::from_rug(tmp_asin),
                expected_acos: RefResult::from_rug(tmp_acos),
            }
        },
        pb,
    );
}

fn gen_args() -> Vec<f32> {
    let mut rng = create_prng();

    let mut args = Vec::new();
    for e in -126..=-1 {
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

    for i in 1..=1000 {
        let x = (i as f32) / 1000.0;
        args.push(x);
        args.push(-x);
    }

    args
}
