use rand::Rng as _;

use super::{calc_error_ulp, mkfloat};
use crate::data::create_prng;

#[test]
fn test_log() {
    let mut max_error: f32 = 0.0;
    test_log_with(|x| {
        let expected = fpmath::log(f64::from(x));
        let actual = fpmath::log(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.9, "log({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max log error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_log2() {
    let mut max_error: f32 = 0.0;
    test_log_with(|x| {
        let expected = fpmath::log2(f64::from(x));
        let actual = fpmath::log2(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.9, "log2({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max log2 error = {max_error}");
    assert!(max_error > 0.5);
}

#[test]
fn test_log10() {
    let mut max_error: f32 = 0.0;
    test_log_with(|x| {
        let expected = fpmath::log10(f64::from(x));
        let actual = fpmath::log10(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.9, "log10({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max log10 error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_log_with(mut f: impl FnMut(f32)) {
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

#[test]
fn test_log_1p() {
    let mut max_error: f32 = 0.0;
    test_log1p_with(|x| {
        let expected = fpmath::log_1p(f64::from(x));
        let actual = fpmath::log_1p(x);

        let err = calc_error_ulp(actual, expected);
        max_error = max_error.max(err);

        assert!(err < 0.9, "log_1p({x:e}) = {actual:e} (error = {err} ULP)");
    });
    eprintln!("max log_1p error = {max_error}");
    assert!(max_error > 0.5);
}

fn test_log1p_with(mut f: impl FnMut(f32)) {
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

    // 1 < x < 0
    for e in -126..=-1 {
        f(mkfloat(0, e, true));
        f(mkfloat(u32::MAX, e, true));

        for _ in 0..1000 {
            let m = rng.random::<u32>();
            f(mkfloat(m, e, true));
        }
    }

    // subnormals
    for i in 0..23 {
        f(f32::from_bits(1 << i));
        f(-f32::from_bits(1 << i));
        f(f32::from_bits((1 << (i + 1)) - 1));
        f(-f32::from_bits((1 << (i + 1)) - 1));
    }
}
