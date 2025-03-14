use rand::Rng as _;

use super::{calc_error_ulp, mkfloat};
use crate::create_prng;

#[test]
fn test_cbrt() {
    let mut max_error: f32 = 0.0;
    test_with(|x| {
        let expected = fpmath::cbrt(f64::from(x));
        let actual = fpmath::cbrt(x);
        assert_eq!(fpmath::cbrt(-x), -actual);

        let err = calc_error_ulp(actual, expected);
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
        f(mkfloat(u32::MAX, e, false));

        for _ in 0..10000 {
            let m = rng.random::<u32>();
            f(mkfloat(m, e, false));
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
        f(f32::from_bits((1 << (i + 1)) - 1));
    }
}
