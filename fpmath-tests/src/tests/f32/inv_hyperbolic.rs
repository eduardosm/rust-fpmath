use rand::Rng as _;

use super::{calc_error_ulp, mkfloat};
use crate::data::create_prng;

#[test]
fn test_asinh() {
    let mut max_error: f32 = 0.0;
    test_asinh_with(|x| {
        let expected = fpmath::asinh(f64::from(x));
        let actual = fpmath::asinh(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.9, "asinh({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max asinh error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_asinh_with(mut f: impl FnMut(f32)) {
    let mut rng = create_prng();

    for e in -126..=127 {
        f(mkfloat(0, e, false));
        f(mkfloat(0, e, true));
        f(mkfloat(u32::MAX, e, false));
        f(mkfloat(u32::MAX, e, true));

        for _ in 0..10000 {
            let m = rng.random::<u32>();
            let s = rng.random::<bool>();
            f(mkfloat(m, e, s));
        }
    }
}

#[test]
fn test_acosh() {
    let mut max_error: f32 = 0.0;
    test_acosh_with(|x| {
        let expected = fpmath::acosh(f64::from(x));
        let actual = fpmath::acosh(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.9, "acosh({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max acosh error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_acosh_with(mut f: impl FnMut(f32)) {
    let mut rng = create_prng();

    for e in 0..=127 {
        f(mkfloat(0, e, false));
        f(mkfloat(u32::MAX, e, false));

        for _ in 0..10000 {
            let m = rng.random::<u32>();
            f(mkfloat(m, e, false));
        }
    }
}

#[test]
fn test_atanh() {
    let mut max_error: f32 = 0.0;
    test_atanh_with(|x| {
        let expected = fpmath::atanh(f64::from(x));
        let actual = fpmath::atanh(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.9, "atanh({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max atanh error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_atanh_with(mut f: impl FnMut(f32)) {
    let mut rng = create_prng();

    for e in -126..=-1 {
        f(mkfloat(0, e, false));
        f(mkfloat(0, e, true));
        f(mkfloat(u32::MAX, e, false));
        f(mkfloat(u32::MAX, e, true));

        for _ in 0..10000 {
            let m = rng.random::<u32>();
            let s = rng.random::<bool>();
            f(mkfloat(m, e, s));
        }
    }
}
