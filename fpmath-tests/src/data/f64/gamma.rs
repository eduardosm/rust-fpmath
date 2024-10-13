use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

pub(crate) fn gen_data_t(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_tgamma",
        gen_args,
        |x| {
            let tmp = rug::Float::with_val(53 * 2, x).gamma();

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
    let pi = rug::Float::with_val(53 * 2, rug::float::Constant::Pi);

    generate_data(
        "f64_lgamma",
        gen_args,
        |x| {
            let mut sign = 1;
            let tmp = if x < 0.5 {
                let bigx = rug::Float::with_val(53 * 2, x);

                let lgamma_omx = rug::Float::with_val(53 * 2, 1u8 - &bigx).ln_gamma();

                let pi_sinpi = &pi / bigx.sin_pi();
                if x < 0.0 && x.fract() == 0.0 {
                    sign = 0;
                } else if pi_sinpi.is_sign_negative() {
                    sign = -1;
                }

                pi_sinpi.abs().ln() - lgamma_omx
            } else {
                rug::Float::with_val(53 * 2, x).ln_gamma()
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
            let m = rng.gen::<u64>();
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
