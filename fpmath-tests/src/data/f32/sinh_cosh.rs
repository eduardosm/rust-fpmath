use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::{create_prng, generate_data};

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f32_sinh_cosh",
        gen_args,
        |x| {
            let mut tmp_arg = dev_mpfr::Mpfr::new();
            tmp_arg.set_prec(24);
            tmp_arg.set_f32(x, dev_mpfr::Rnd::N);

            let mut tmp_sin = dev_mpfr::Mpfr::new();
            tmp_sin.set_prec(24 * 2);

            let mut tmp_cos = dev_mpfr::Mpfr::new();
            tmp_cos.set_prec(24 * 2);

            tmp_arg.sinh_cosh(&mut tmp_sin, &mut tmp_cos, dev_mpfr::Rnd::N);

            super::SinCosData {
                x,
                expected_sin: RefResult::from_mpfr(&mut tmp_sin),
                expected_cos: RefResult::from_mpfr(&mut tmp_cos),
            }
        },
        pb,
    );
}

pub(super) fn gen_args() -> Vec<f32> {
    let mut rng = create_prng();

    let mut args = Vec::new();
    for e in -126..=9 {
        args.push(mkfloat(0, e, false));
        args.push(mkfloat(0, e, true));
        args.push(mkfloat(u32::MAX, e, false));
        args.push(mkfloat(u32::MAX, e, true));

        for _ in 0..1000 {
            let m = rng.gen::<u32>();
            let s = rng.gen::<bool>();
            args.push(mkfloat(m, e, s));
        }
    }

    for arg in -1000..=1000 {
        args.push(arg as f32);
    }

    args
}
