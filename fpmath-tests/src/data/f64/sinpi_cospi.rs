use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

pub(super) fn gen_args() -> Vec<f64> {
    let mut rng = create_prng();

    let mut args = Vec::new();
    for e in -1022..=1023 {
        args.push(mkfloat(0, e, false));
        args.push(mkfloat(0, e, true));
        args.push(mkfloat(u64::MAX, e, false));
        args.push(mkfloat(u64::MAX, e, true));

        for _ in 0..1000 {
            let m = rng.gen::<u64>();
            let s = rng.gen::<bool>();
            args.push(mkfloat(m, e, s));
        }
    }

    for arg in -20000..=20000 {
        args.push((arg as f64) * 0.5);
    }

    args
}

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    let mut k_2 = dev_mpfr::Mpfr::new();
    k_2.set_prec(128);
    k_2.set_ui(2, dev_mpfr::Rnd::N);

    let mut conv = dev_mpfr::Mpfr::new();
    conv.set_prec(2048);
    conv.const_pi(dev_mpfr::Rnd::N);

    generate_data(
        "f64_sinpi_cospi",
        gen_args,
        |x| {
            let mut tmp_arg = dev_mpfr::Mpfr::new();
            tmp_arg.set_prec(2048);
            tmp_arg.set_f64(x, dev_mpfr::Rnd::N);
            tmp_arg.fmod(None, Some(&k_2), dev_mpfr::Rnd::N);
            tmp_arg.mul(Some(&conv), None, dev_mpfr::Rnd::N);

            let mut tmp_sin = dev_mpfr::Mpfr::new();
            tmp_sin.set_prec(53 * 2);

            let mut tmp_cos = dev_mpfr::Mpfr::new();
            tmp_cos.set_prec(53 * 2);

            tmp_arg.sin_cos(&mut tmp_sin, &mut tmp_cos, dev_mpfr::Rnd::N);

            super::SinCosData {
                x,
                expected_sin: RefResult::from_mpfr(&mut tmp_sin),
                expected_cos: RefResult::from_mpfr(&mut tmp_cos),
            }
        },
        pb,
    );
}
