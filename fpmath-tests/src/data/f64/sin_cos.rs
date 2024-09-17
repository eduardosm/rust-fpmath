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

    for arg in -1000..=1000 {
        args.push(arg as f64);
    }

    // Problematic value in
    // "ARGUMENT REDUCTION FOR HUGE ARGUMENTS: Good to the Last Bit"
    args.push(1.0e22);

    // Magic values from rust-libm tests
    args.push(3.141592025756836);
    args.push(3.141592033207416);
    args.push(3.141592144966125);
    args.push(3.141592979431152);
    args.push(-3054214.5490637687);
    args.push(917340800458.2274);

    let f2s = [1.0, -1.0, (1 << 20) as f64, -((1 << 20) as f64)];
    let f3s = [
        1.0, 1.1, 1.01, 1.001, 1.0001, 1.00001, 1.000001, 1.0000001, 0.9, 0.99, 0.999, 0.9999,
        0.99999, 0.999999, 0.9999999,
    ];

    for f1 in 1..=100 {
        for f2 in f2s {
            for f3 in f3s {
                args.push(std::f64::consts::FRAC_PI_8 * (f1 as f64) * f2 * f3);
            }
        }
    }

    args
}

pub(crate) fn gen_data(pb: indicatif::ProgressBar) {
    generate_data(
        "f64_sin_cos",
        gen_args,
        |x| {
            let mut tmp_arg = dev_mpfr::Mpfr::new();
            tmp_arg.set_prec(53);
            tmp_arg.set_f64(x, dev_mpfr::Rnd::N);

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
