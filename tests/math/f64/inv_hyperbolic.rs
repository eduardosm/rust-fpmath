use rand::Rng as _;

use super::{calc_error_ulp, mkfloat, RUG_PREC};
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
            let m = rng.random::<u64>();
            let s = rng.random::<bool>();
            f(mkfloat(m, e, s));
        }
    }
    for e in -1022..=1023 {
        f(mkfloat(0, e, false));
        f(mkfloat(0, e, true));
        f(mkfloat(u64::MAX, e, false));
        f(mkfloat(u64::MAX, e, true));

        for _ in 0..1000 {
            let m = rng.random::<u64>();
            let s = rng.random::<bool>();
            f(mkfloat(m, e, s));
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
        f(mkfloat(0, e, false));
        f(mkfloat(u64::MAX, e, false));

        for _ in 0..9000 {
            let m = rng.random::<u64>();
            f(mkfloat(m, e, false));
        }
    }

    for e in 0..=1023 {
        f(mkfloat(0, e, false));
        f(mkfloat(u64::MAX, e, false));

        for _ in 0..1000 {
            let m = rng.random::<u64>();
            f(mkfloat(m, e, false));
        }
    }
}

#[test]
fn test_atanh() {
    let mut max_error: f64 = 0.0;
    test_atanh_with(|x| {
        let expected = rug::Float::with_val(RUG_PREC, x).atanh();
        let actual = fpmath::atanh(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.9, "atanh({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max atanh error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_atanh_with(mut f: impl FnMut(f64)) {
    let mut rng = create_prng();

    for e in -1022..=-1 {
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
}
