use rand::Rng as _;

use super::{calc_error_ulp, mkfloat, RUG_PREC};
use crate::create_prng;

#[test]
fn test_cbrt() {
    let mut max_error: f64 = 0.0;
    test_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).cbrt();
        let actual = fpmath::cbrt(x);
        assert_eq!(fpmath::cbrt(-x), -actual);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.9, "cbrt({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max cbrt error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_with(mut f: impl FnMut(f64)) {
    let mut rng = create_prng();

    for e in -1022..=1023 {
        f(mkfloat(0, e, false));
        f(mkfloat(u64::MAX, e, false));

        for _ in 0..5000 {
            let m = rng.random::<u64>();
            f(mkfloat(m, e, false));
        }
    }
    for e in -1022..=1023 {
        for _ in 0..5000 {
            let m = rng.random::<u64>();
            f(mkfloat(m, e, false));
        }
    }

    for arg in 1..=10000 {
        f(arg as f64);
    }

    f(f64::MIN_POSITIVE);
    f(f64::MAX);

    // subnormals
    for i in 0..52 {
        f(f64::from_bits(1 << i));
        f(f64::from_bits((1 << (i + 1)) - 1));
    }
}
