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
    generate_data(
        "f32_sinpi_cospi",
        gen_args,
        |x| {
            let tmp_arg = rug::Float::with_val(24, x);
            let tmp_sin = rug::Float::with_val(24 * 2, tmp_arg.sin_pi_ref());
            let tmp_cos = rug::Float::with_val(24 * 2, tmp_arg.cos_pi_ref());

            super::SinCosData {
                x,
                expected_sin: RefResult::from_rug(tmp_sin),
                expected_cos: RefResult::from_rug(tmp_cos),
            }
        },
        pb,
    );
}
