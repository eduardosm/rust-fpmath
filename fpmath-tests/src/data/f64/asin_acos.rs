use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

#[derive(bincode::Encode, bincode::Decode)]
pub(crate) struct Data {
    pub(crate) x: f64,
    pub(crate) expected_asin: RefResult,
    pub(crate) expected_acos: RefResult,
}

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_asin_acos",
        gen_args,
        |x| {
            let bigx = rug::Float::with_val(53, x);
            let tmp_asin = rug::Float::with_val(53 * 2, bigx.asin_ref());
            let tmp_acos = rug::Float::with_val(53 * 2, bigx.acos_ref());

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
    generate_data(
        "f64_asind_acosd",
        gen_args,
        |x| {
            let bigx = rug::Float::with_val(53, x);
            let tmp_asin = rug::Float::with_val(53 * 2, bigx.asin_u_ref(360));
            let tmp_acos = rug::Float::with_val(53 * 2, bigx.acos_u_ref(360));

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
        "f64_asinpi_acospi",
        gen_args,
        |x| {
            let bigx = rug::Float::with_val(53, x);
            let tmp_asin = rug::Float::with_val(53 * 2, bigx.asin_pi_ref());
            let tmp_acos = rug::Float::with_val(53 * 2, bigx.acos_pi_ref());

            Data {
                x,
                expected_asin: RefResult::from_rug(tmp_asin),
                expected_acos: RefResult::from_rug(tmp_acos),
            }
        },
        pb,
    );
}

fn gen_args() -> Vec<f64> {
    let mut rng = create_prng();

    let mut args = Vec::new();
    for e in -1022..=-1 {
        args.push(mkfloat(0, e, false));
        args.push(mkfloat(0, e, true));
        args.push(mkfloat(u64::MAX, e, false));
        args.push(mkfloat(u64::MAX, e, true));

        for _ in 0..10000 {
            let m = rng.random::<u64>();
            let s = rng.random::<bool>();
            args.push(mkfloat(m, e, s));
        }
    }

    for i in 1..=1000 {
        let x = (i as f64) / 1000.0;
        args.push(x);
        args.push(-x);
    }

    args
}
