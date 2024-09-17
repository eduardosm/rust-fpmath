use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

pub(super) fn gen_args() -> Vec<f32> {
    let mut rng = create_prng();

    let mut args = Vec::new();
    for e in -126..=127 {
        args.push(mkfloat(0, e, false));
        args.push(mkfloat(0, e, true));
        args.push(mkfloat(u32::MAX, e, false));
        args.push(mkfloat(u32::MAX, e, true));

        for _ in 0..5000 {
            let m = rng.gen::<u32>();
            let s = rng.gen::<bool>();
            args.push(mkfloat(m, e, s));
        }
    }

    for arg in -20000..=20000 {
        args.push((arg as f32) * 0.5);
    }

    args
}

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    let mut k_360 = dev_mpfr::Mpfr::new();
    k_360.set_prec(64);
    k_360.set_ui(360, dev_mpfr::Rnd::N);

    let mut conv = dev_mpfr::Mpfr::new();
    conv.set_prec(512);
    conv.const_pi(dev_mpfr::Rnd::N);
    conv.div_f64(None, 180.0, dev_mpfr::Rnd::N);

    generate_data(
        "f32_sind_cosd",
        gen_args,
        |x| {
            let mut tmp_arg = dev_mpfr::Mpfr::new();
            tmp_arg.set_prec(512);
            tmp_arg.set_f32(x, dev_mpfr::Rnd::N);
            tmp_arg.fmod(None, Some(&k_360), dev_mpfr::Rnd::N);
            tmp_arg.mul(Some(&conv), None, dev_mpfr::Rnd::N);

            let mut tmp_sin = dev_mpfr::Mpfr::new();
            tmp_sin.set_prec(24 * 2);

            let mut tmp_cos = dev_mpfr::Mpfr::new();
            tmp_cos.set_prec(24 * 2);

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
