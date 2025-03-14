use rand::Rng as _;

use super::{calc_error_ulp, mkfloat, purify, purify2, RUG_PREC};
use crate::create_prng;

#[test]
fn test_sinh_cosh() {
    let mut max_sin1_error: f64 = 0.0;
    let mut max_sin2_error: f64 = 0.0;
    let mut max_cos1_error: f64 = 0.0;
    let mut max_cos2_error: f64 = 0.0;
    test_with(|x| {
        let (expected_sin, expected_cos) =
            rug::Float::with_val(RUG_PREC, x).sinh_cosh(rug::Float::new(RUG_PREC));

        let actual_sin1 = fpmath::sinh(x);
        let actual_cos1 = fpmath::cosh(x);
        let (actual_sin2, actual_cos2) = fpmath::sinh_cosh(x);
        assert_eq!(purify(fpmath::sinh(-x)), purify(-actual_sin1));
        assert_eq!(purify(fpmath::cosh(-x)), purify(actual_cos1));
        assert_eq!(
            purify2(fpmath::sinh_cosh(-x)),
            purify2((-actual_sin2, actual_cos2))
        );

        let sin1_err = calc_error_ulp(actual_sin1, expected_sin.clone());
        let sin2_err = calc_error_ulp(actual_sin2, expected_sin);
        let cos1_err = calc_error_ulp(actual_cos1, expected_cos.clone());
        let cos2_err = calc_error_ulp(actual_cos2, expected_cos);

        max_sin1_error = max_sin1_error.max(sin1_err);
        max_sin2_error = max_sin2_error.max(sin2_err);
        max_cos1_error = max_cos1_error.max(cos1_err);
        max_cos2_error = max_cos2_error.max(cos2_err);

        assert!(
            sin1_err < 0.9,
            "sinh({x:e}) = {actual_sin1:e} (error = {sin1_err} ULP)",
        );
        assert!(
            sin2_err < 0.9,
            "sinh({x:e}) = {actual_sin2:e} (error = {sin2_err} ULP)",
        );

        assert!(
            cos1_err < 0.9,
            "cosh({x:e}) = {actual_cos1:e} (error = {cos1_err} ULP)",
        );
        assert!(
            cos2_err < 0.9,
            "cosh({x:e}) = {actual_cos2:e} (error = {cos2_err} ULP)",
        );
    });
    eprintln!("max sinh1 error = {max_sin1_error}");
    eprintln!("max sinh2 error = {max_sin2_error}");
    eprintln!("max cosh1 error = {max_cos1_error}");
    eprintln!("max cosh2 error = {max_cos2_error}");
    assert!(max_sin1_error > 0.5);
    assert!(max_sin2_error > 0.5);
    assert!(max_cos1_error > 0.5);
    assert!(max_cos2_error > 0.5);
}

#[test]
fn test_tanh() {
    let mut max_error: f64 = 0.0;
    test_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).tanh();
        let actual = fpmath::tanh(x);
        assert_eq!(fpmath::tanh(-x), -actual);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.9, "tanh({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max tanh error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_with(mut f: impl FnMut(f64)) {
    let mut rng = create_prng();

    for e in -1022..-200 {
        f(mkfloat(0, e, false));
        f(mkfloat(u64::MAX, e, false));

        for _ in 0..100 {
            let m = rng.random::<u64>();
            f(mkfloat(m, e, false));
        }
    }
    for e in -200..=12 {
        f(mkfloat(0, e, false));
        f(mkfloat(u64::MAX, e, false));

        for _ in 0..1000 {
            let m = rng.random::<u64>();
            f(mkfloat(m, e, false));
        }
    }

    for arg in 1..=1000 {
        f(arg as f64);
    }
}
