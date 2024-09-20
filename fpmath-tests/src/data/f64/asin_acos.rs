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
            let mut tmp_arg = dev_mpfr::Mpfr::new();
            tmp_arg.set_prec(53);
            tmp_arg.set_f64(x, dev_mpfr::Rnd::N);

            let mut tmp_asin = dev_mpfr::Mpfr::new();
            tmp_asin.set_prec(53 * 2);

            let mut tmp_acos = dev_mpfr::Mpfr::new();
            tmp_acos.set_prec(53 * 2);

            tmp_asin.asin(Some(&tmp_arg), dev_mpfr::Rnd::N);
            tmp_acos.acos(Some(&tmp_arg), dev_mpfr::Rnd::N);

            Data {
                x,
                expected_asin: RefResult::from_mpfr(&mut tmp_asin),
                expected_acos: RefResult::from_mpfr(&mut tmp_acos),
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
        "f64_asind_acosd",
        gen_args,
        |x| {
            let mut tmp_arg = dev_mpfr::Mpfr::new();
            tmp_arg.set_prec(53);
            tmp_arg.set_f64(x, dev_mpfr::Rnd::N);

            let mut tmp_asin = dev_mpfr::Mpfr::new();
            tmp_asin.set_prec(53 * 3);

            let mut tmp_acos = dev_mpfr::Mpfr::new();
            tmp_acos.set_prec(53 * 3);

            tmp_asin.asin(Some(&tmp_arg), dev_mpfr::Rnd::N);
            tmp_acos.acos(Some(&tmp_arg), dev_mpfr::Rnd::N);

            tmp_asin.mul(None, Some(&conv), dev_mpfr::Rnd::N);
            tmp_acos.mul(None, Some(&conv), dev_mpfr::Rnd::N);

            Data {
                x,
                expected_asin: RefResult::from_mpfr(&mut tmp_asin),
                expected_acos: RefResult::from_mpfr(&mut tmp_acos),
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
        "f64_asinpi_acospi",
        gen_args,
        |x| {
            let mut tmp_arg = dev_mpfr::Mpfr::new();
            tmp_arg.set_prec(53);
            tmp_arg.set_f64(x, dev_mpfr::Rnd::N);

            let mut tmp_asin = dev_mpfr::Mpfr::new();
            tmp_asin.set_prec(53 * 3);

            let mut tmp_acos = dev_mpfr::Mpfr::new();
            tmp_acos.set_prec(53 * 3);

            tmp_asin.asin(Some(&tmp_arg), dev_mpfr::Rnd::N);
            tmp_acos.acos(Some(&tmp_arg), dev_mpfr::Rnd::N);

            tmp_asin.div(None, Some(&conv), dev_mpfr::Rnd::N);
            tmp_acos.div(None, Some(&conv), dev_mpfr::Rnd::N);

            Data {
                x,
                expected_asin: RefResult::from_mpfr(&mut tmp_asin),
                expected_acos: RefResult::from_mpfr(&mut tmp_acos),
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
            let m = rng.gen::<u64>();
            let s = rng.gen::<bool>();
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
