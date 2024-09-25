use rand::Rng as _;

use super::{mkfloat, RefResult};
use crate::data::create_prng;

#[test]
fn test_cbrt() {
    let mut max_error: f32 = 0.0;
    test_with(|x| {
        let expected = RefResult::from_f64(fpmath::cbrt(f64::from(x)));
        let actual = fpmath::cbrt(x);

        let err = expected.calc_error(actual);
        max_error = max_error.max(err);

        assert!(err < 0.9, "cbrt({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max cbrt error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_with(mut f: impl FnMut(f32)) {
    let mut rng = create_prng();

    for e in -126..=127 {
        f(mkfloat(0, e, false));
        f(mkfloat(0, e, true));
        f(mkfloat(u32::MAX, e, false));
        f(mkfloat(u32::MAX, e, true));

        for _ in 0..10000 {
            let m = rng.gen::<u32>();
            let s = rng.gen::<bool>();
            f(mkfloat(m, e, s));
        }
    }

    for arg in 1..=10000 {
        f(arg as f32);
    }

    f(f32::MIN_POSITIVE);
    f(f32::MAX);

    // subnormals
    for i in 0..23 {
        f(f32::from_bits(1 << i));
        f(-f32::from_bits(1 << i));
        f(f32::from_bits((1 << (i + 1)) - 1));
        f(-f32::from_bits((1 << (i + 1)) - 1));
    }
}
