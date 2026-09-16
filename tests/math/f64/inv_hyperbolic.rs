use rand::RngExt as _;

use super::{RUG_PREC, calc_error_ulp, mk_normal, purify};
use crate::create_prng;

#[test]
fn test_asinh() {
    let mut max_error: f64 = 0.0;
    test_asinh_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).asinh();
        let actual = fpmath::asinh(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.9, "asinh({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max asinh error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_asinh_with(mut f: impl FnMut(f64)) {
    let mut rng = create_prng();

    for e in -100..=100 {
        for _ in 0..9000 {
            let m = super::gen_mantissa(&mut rng);
            let s = rng.random::<bool>();
            f(mk_normal(m, e, s));
        }
    }
    for e in -1022..=1023 {
        f(mk_normal(0, e, false));
        f(mk_normal(0, e, true));
        f(mk_normal(super::MAX_MANTISSA, e, false));
        f(mk_normal(super::MAX_MANTISSA, e, true));

        for _ in 0..1000 {
            let m = super::gen_mantissa(&mut rng);
            let s = rng.random::<bool>();
            f(mk_normal(m, e, s));
        }
    }
}

#[test]
fn test_acosh() {
    let mut max_error: f64 = 0.0;
    test_acosh_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).acosh();
        let actual = fpmath::acosh(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.9, "acosh({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max acosh error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_acosh_with(mut f: impl FnMut(f64)) {
    let mut rng = create_prng();

    for e in 0..=100 {
        f(mk_normal(0, e, false));
        f(mk_normal(super::MAX_MANTISSA, e, false));

        for _ in 0..9000 {
            let m = super::gen_mantissa(&mut rng);
            f(mk_normal(m, e, false));
        }
    }

    for e in 0..=1023 {
        f(mk_normal(0, e, false));
        f(mk_normal(super::MAX_MANTISSA, e, false));

        for _ in 0..1000 {
            let m = super::gen_mantissa(&mut rng);
            f(mk_normal(m, e, false));
        }
    }
}

#[test]
fn test_atanh() {
    let mut max_error: f64 = 0.0;
    test_atanh_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).atanh();
        let actual = fpmath::atanh(x);
        assert_eq!(purify(fpmath::atanh(-x)), purify(-actual));

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.9, "atanh({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max atanh error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_atanh_with(mut f: impl FnMut(f64)) {
    let mut rng = create_prng();

    // Near ±1
    for m0 in 0..(1 << 23) {
        f(mk_normal(super::MAX_MANTISSA - m0, -1, false));
    }

    for e in -1022..=-1 {
        f(mk_normal(0, e, false));
        f(mk_normal(super::MAX_MANTISSA, e, false));

        for _ in 0..10000 {
            let m = super::gen_mantissa(&mut rng);
            let s = rng.random::<bool>();
            f(mk_normal(m, e, s));
        }
    }
}
