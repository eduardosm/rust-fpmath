use rand::Rng as _;

use super::{mkfloat, RefResult, RUG_PREC};
use crate::data::{create_prng, generate_data};

pub(crate) fn gen_data_t(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_tgamma",
        gen_args,
        |x| {
            let tmp = rug::Float::with_val(RUG_PREC, x).gamma();

            super::OneArgData {
                x,
                expected: RefResult::from_rug(tmp),
            }
        },
        pb,
    );
}

#[derive(bincode::Encode, bincode::Decode)]
pub(crate) struct LgammaData {
    pub(crate) x: f64,
    pub(crate) expected: RefResult,
    pub(crate) expected_sign: i8,
}

pub(crate) fn gen_data_l(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_lgamma",
        gen_args,
        |x| {
            let (tmp, ord) = rug::Float::with_val(RUG_PREC, x).ln_abs_gamma();
            let sign = if x < 0.0 && x.fract() == 0.0 {
                0
            } else {
                ord as i8
            };

            LgammaData {
                x,
                expected: RefResult::from_rug(tmp),
                expected_sign: sign,
            }
        },
        pb,
    );
}

fn gen_args() -> Vec<f64> {
    let mut rng = create_prng();

    let mut args = Vec::new();
    for e in -1022..=1023 {
        for _ in 0..3000 {
            let m = rng.random::<u64>();
            args.push(mkfloat(m, e, false));
            args.push(mkfloat(m, e, true));
        }
    }

    for i in 0..20000 {
        let x = (i as f64) / 100.0;
        args.push(x);
        args.push(-x);
    }

    args
}
