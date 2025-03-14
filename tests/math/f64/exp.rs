use rand::Rng as _;

use super::{calc_error_ulp, mkfloat, select_threshold, RUG_PREC};
use crate::create_prng;

#[test]
fn test_exp() {
    let mut max_exp_error: f64 = 0.0;
    let mut max_expm1_error: f64 = 0.0;
    test_with(|x| {
        let expected_exp = rug::Float::with_val(RUG_PREC, x).exp();
        let expected_expm1 = rug::Float::with_val(RUG_PREC, x).exp_m1();

        let actual_exp = fpmath::exp(x);
        let actual_expm1 = fpmath::exp_m1(x);

        let exp_err = calc_error_ulp(actual_exp, expected_exp);
        let expm1_err = calc_error_ulp(actual_expm1, expected_expm1);

        max_exp_error = max_exp_error.max(exp_err);
        max_expm1_error = max_expm1_error.max(expm1_err);

        let exp_threshold = select_threshold(actual_exp, 0.9, 1.9);
        assert!(
            exp_err < exp_threshold,
            "exp({x:e}) = {actual_exp:e} (error = {exp_err} ULP)",
        );

        let expm1_threshold = select_threshold(actual_expm1, 0.9, 1.9);
        assert!(
            expm1_err < expm1_threshold,
            "expm1({x:e}) = {actual_expm1:e} (error = {expm1_err} ULP)",
        );
    });
    eprintln!("max exp error = {max_exp_error}");
    eprintln!("max expm1 error = {max_expm1_error}");
    assert!(max_exp_error > 0.5);
    assert!(max_expm1_error > 0.5);
}

#[test]
fn test_exp2() {
    let mut max_error: f64 = 0.0;
    test_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).exp2();
        let actual = fpmath::exp2(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.9, 1.9);
        assert!(
            err < threshold,
            "exp2({x:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max exp2 error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_exp10() {
    let mut max_error: f64 = 0.0;
    test_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).exp10();
        let actual = fpmath::exp10(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        let threshold = select_threshold(actual, 0.9, 1.9);
        assert!(
            err < threshold,
            "exp10({x:e}) = {actual:e} (error = {err} ULP)",
        );
    });
    eprintln!("max exp10 error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_with(mut f: impl FnMut(f64)) {
    let mut rng = create_prng();

    for e in -1022..=12 {
        f(mkfloat(0, e, false));
        f(mkfloat(0, e, true));
        f(mkfloat(u64::MAX, e, false));
        f(mkfloat(u64::MAX, e, true));

        for _ in 0..10000 {
            let m = rng.random::<u64>();
            let s = rng.random::<bool>();
            f(mkfloat(m, e, s));
        }
    }

    for arg in -1100..=1100 {
        f(arg as f64);
    }
}
