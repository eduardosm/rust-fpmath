use super::{RUG_PREC, calc_error_ulp, mk_normal};
use crate::create_prng;

#[test]
fn test_exp() {
    let mut max_error: f64 = 0.0;
    test_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).exp();
        let actual = fpmath::exp(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.51, "exp({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_exp_m1() {
    let mut max_error: f64 = 0.0;
    test_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).exp_m1();
        let actual = fpmath::exp_m1(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.51, "exp_m1({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_exp2() {
    let mut max_error: f64 = 0.0;
    test_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).exp2();
        let actual = fpmath::exp2(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.51, "exp2({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max error = {max_error}");
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

        assert!(err < 0.51, "exp10({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_with(mut f: impl FnMut(f64)) {
    let mut rng = create_prng();

    for e in -1022..=12 {
        f(mk_normal(0, e, false));
        f(mk_normal(0, e, true));
        f(mk_normal(super::MAX_MANTISSA, e, false));
        f(mk_normal(super::MAX_MANTISSA, e, true));

        for _ in 0..50000 {
            let m = super::gen_mantissa(&mut rng);
            f(mk_normal(m, e, false));
            f(mk_normal(m, e, true));
        }
    }

    for arg in -1100..=1100 {
        f(arg as f64);
    }
}
